/*
 * Copyright (c) 2014, 2015 Manu T S
 *
 * This is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this software.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "framing.h"
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <volk/volk.h>
#include <volk/volk_malloc.h>
#include <string.h>

#define ADD_NULL_CARRIERS   false

gr_complex CHANNEL_I[2][2] =
{
  {gr_complex(1.0, 0.0),
   gr_complex(0.0, 0.0)},
  {gr_complex(0.0, 0.0),
   gr_complex(1.0, 0.0)}
};

gr_complex BPSK_CONSTELLATION[] = 
{
  gr_complex(-1.0, 0.0),
  gr_complex(1.0,  0.0),
};
#define BPSK_CONSTELLATION_SIZE 2
gr_complex QPSK_CONSTELLATION[] = 
{
  gr_complex(0.5, 0.5),
  gr_complex(-0.5, 0.5),
  gr_complex(-0.5, -0.5),
  gr_complex(0.5, -0.5)
};
#define QPSK_CONSTELLATION_SIZE 4
#define CONSTELLATION BPSK_CONSTELLATION
#define CONSTELLATION_SIZE BPSK_CONSTELLATION_SIZE

inline gr_complex
liquid_cexpjf(float theta)
{
  return std::polar(1.0f, theta);
}

inline float
cabsf(gr_complex z)
{
  return std::abs(z);
}

inline float
cargf(gr_complex z)
{
  return std::arg(z);
}

inline float
fabsf(float x)
{
  return std::abs(x);
}

inline gr_complex
conjf(gr_complex z)
{
  return std::conj(z);
}

framegen::framegen(unsigned int _M,
                   unsigned int _cp_len,
                   unsigned int _num_streams,
                   unsigned int _num_access_codes,
                   unsigned char * const &_p,
                   msequence const &_ms_S0,
                   std::vector<msequence> const &_ms_S1) 
{
  // assign values for variables
  M = _M;
  cp_len = _cp_len;
  symbol_len = M + cp_len;
  num_streams = _num_streams;
  num_access_codes = _num_access_codes;

  // resize vectors to num stream
  W = (gr_complex ***) malloc 
      (sizeof(gr_complex [M][num_streams][num_streams]));
  // initialize the precoding to unity
  gr_complex I(1.0, 0.0);
  std::fill((gr_complex *)W, ((gr_complex *)W) + M*num_streams*num_streams, I);
  ifft.resize(num_streams);
  X.resize(num_streams);
  x.resize(num_streams);
  S1.resize(num_streams);
  s1.resize(num_streams);
  p  = (unsigned char *) malloc (sizeof(unsigned char)*M);
  memmove(p, _p, sizeof(unsigned char)*M);
  S0 = (gr_complex *) malloc (sizeof(gr_complex)*M);
  s0 = (gr_complex *) malloc (sizeof(gr_complex)*M);
  // sanity check for subcarrier allocation
  ofdmframe_validate_sctype(p,
                            M,
                            &M_null,
                            &M_pilot,
                            &M_data);
  // initialize S0 sequence
  ofdmframe_init_S0(p,
                    M,
                    S0,
                    s0,
                    _ms_S0);
  g_data = 1.0f / sqrtf(M_pilot + M_data);

  // get volk alignment for volk_malloc
  size_t volk_alignment = volk_get_alignment();
  for(unsigned int i = 0; i < num_streams; i++)
  {
    // memory allocation
    // NOTE: all to be freed in destructor
    S1[i] = (gr_complex *) malloc 
            (sizeof(gr_complex)*
             M*num_access_codes);
    s1[i] = (gr_complex *) malloc 
            (sizeof(gr_complex)*
             M*num_access_codes);
    X[i] = (std::complex<float> *) fftwf_malloc
           (sizeof(std::complex<float>)*M);
    x[i] = (std::complex<float> *) fftwf_malloc
           (sizeof(std::complex<float>)*M);
    assert(volk_is_aligned(X[i]));
    assert(volk_is_aligned(x[i]));
    ifft[i] = fftwf_plan_dft_1d(M,
				reinterpret_cast<fftwf_complex *>(X[i]),
				reinterpret_cast<fftwf_complex *>(x[i]),
                                FFTW_BACKWARD,
                                FFTW_ESTIMATE);
    // initialize S1 sequence
    ofdmframe_init_S1(p,
                      M,
                      num_access_codes,
                      S1[i],
                      s1[i],
                      _ms_S1[i]);
  }
  volk_buff_fc1 = (std::complex<float> *) volk_malloc 
                  (sizeof(std::complex<float>)*M,
                   volk_alignment);
  volk_buff_fc2 = (std::complex<float> *) volk_malloc 
                  (sizeof(std::complex<float>)*M,
                   volk_alignment);
  volk_buff_fc3 = (std::complex<float> *) volk_malloc 
                  (sizeof(std::complex<float>)*M,
                   volk_alignment);
  volk_buff_f1 = (float *) volk_malloc 
                 (sizeof(float)*
                  (M + cp_len)*3*num_streams,
                  volk_alignment);
}

unsigned int
framegen::get_num_streams()
{
  return num_streams;
}

void
framegen::set_W(gr_complex *** const _W)
{
  // TODO sanity checks on precoding vector
  memmove(W, _W, sizeof(gr_complex[M][num_streams][num_streams]));
}

unsigned int
framegen::write_sync_words(std::vector<std::complex<float> *> tx_buff)
{
  assert(num_streams == tx_buff.size());
  unsigned int sample_index = 0;
  unsigned int total_count = (num_access_codes*num_streams + 1)*
                             (M + cp_len);
  // set tx_buff to 0;
  for(unsigned int stream = 0; stream < num_streams; stream++) {
    std::fill(tx_buff[stream],
              tx_buff[stream] + total_count,
              gr_complex(0.0, 0.0));
  }
  // write S0 onto ch0 and then S1 in TDMA
  memmove(tx_buff[0] + sample_index,
          s0 + M - cp_len,
          sizeof(std::complex<float>)*cp_len);
  sample_index += cp_len;
  memmove(tx_buff[0] + sample_index,
          s0,
          sizeof(std::complex<float>)*M);
  sample_index += M;
  for(unsigned int ac_id = 0; ac_id < num_access_codes; ac_id++) {
    for(unsigned int stream = 0; stream < num_streams; stream++) {
      // add cyclic prefix
      memmove(tx_buff[stream] + sample_index,
              s1[stream] + M*ac_id + M - cp_len,
              sizeof(std::complex<float>)*cp_len);
      sample_index += cp_len;
      // write ofdm symbol
      memmove(tx_buff[stream] + sample_index,
              s1[stream] + M*ac_id,
              sizeof(std::complex<float>)*M);
      sample_index += M;
    }
  }
  // total count of sync samples
  assert(sample_index == total_count);
  return sample_index;
}

void framegen::compute_W(gr_complex *** G)
{
}

unsigned int
framegen::write_mimo_packet(std::vector<gr_complex *> tx_buff)
{
  return symbol_len;
}

framegen::~framegen()
{
  for(unsigned int i = 0; i < num_streams; i++)
  {
    free(S1[i]);
    free(s1[i]);
    fftwf_free(X[i]);
    fftwf_free(x[i]);
    fftwf_destroy_plan(ifft[i]);
  }
  free(p);
  free(W);
  free(S0);
  free(s0);
  volk_free(volk_buff_fc1);
  volk_free(volk_buff_fc2);
  volk_free(volk_buff_fc3);
  volk_free(volk_buff_f1);
}

void framegen::print()
{
    printf("ofdmframegen:\n");
    printf("    num subcarriers     :   %-u\n", M);
    printf("      - NULL            :   %-u\n", M_null);
    printf("      - pilot           :   %-u\n", M_pilot);
    printf("      - data            :   %-u\n", M_data);
    printf("    cyclic prefix len   :   %-u\n", cp_len);
    printf("    ");
    ofdmframe_print_sctype(p, M);
}

framesync::framesync(unsigned int _M,
                     unsigned int _cp_len,
                     unsigned int _num_streams,
                     unsigned int _num_access_codes,
                     unsigned char * const &_p,
                     msequence const &_ms_S0,
                     std::vector<msequence> const &_ms_S1,
		     mimo_callback _callback) 
{
  // assign values for variables
  M = _M;
  cp_len = _cp_len;
  symbol_len = M + cp_len;
  num_streams = _num_streams;
  num_access_codes = _num_access_codes;

  // resize vectors to num stream
  G = (gr_complex ***) malloc 
      (sizeof(gr_complex [M][num_streams][num_streams]));

  fft.resize(num_streams);
  X.resize(num_streams);
  x.resize(num_streams);
  S1.resize(num_streams);
  s1.resize(num_streams);
  p  = (unsigned char *) malloc (sizeof(unsigned char)*M);
  memmove(p, _p, sizeof(unsigned char)*M);
  S0 = (gr_complex *) malloc (sizeof(gr_complex)*M);
  s0 = (gr_complex *) malloc (sizeof(gr_complex)*M);
  // sanity check for subcarrier allocation
  ofdmframe_validate_sctype(p,
                            M,
                            &M_null,
                            &M_pilot,
                            &M_data);
  // initialize S0 sequence
  ofdmframe_init_S0(p,
                    M,
                    S0,
                    s0,
                    _ms_S0);

  // get volk alignment for volk_malloc
  size_t volk_alignment = volk_get_alignment();
  for(unsigned int i = 0; i < num_streams; i++)
  {
    // memory allocation
    // NOTE: all to be freed in destructor
    S1[i] = (gr_complex *) malloc 
            (sizeof(gr_complex)*
             M*num_access_codes);
    s1[i] = (gr_complex *) malloc 
            (sizeof(gr_complex)*
             M*num_access_codes);
    X[i] = (std::complex<float> *) fftwf_malloc
           (sizeof(std::complex<float>)*M);
    x[i] = (std::complex<float> *) fftwf_malloc
           (sizeof(std::complex<float>)*M);
    assert(volk_is_aligned(X[i]));
    assert(volk_is_aligned(x[i]));
    fft[i] = fftwf_plan_dft_1d(M,
			       reinterpret_cast<fftwf_complex *>(x[i]),
			       reinterpret_cast<fftwf_complex *>(X[i]),
                               FFTW_FORWARD,
                               FFTW_ESTIMATE);
    // initialize S1 sequence
    ofdmframe_init_S1(p,
                      M,
                      num_access_codes,
                      S1[i],
                      s1[i],
                      _ms_S1[i]);
  }
  volk_buff_fc1 = (std::complex<float> *) volk_malloc 
                  (sizeof(std::complex<float>)*M,
                   volk_alignment);
  volk_buff_fc2 = (std::complex<float> *) volk_malloc 
                  (sizeof(std::complex<float>)*M,
                   volk_alignment);
  volk_buff_fc3 = (std::complex<float> *) volk_malloc 
                  (sizeof(std::complex<float>)*M,
                   volk_alignment);
  volk_buff_f1 = (float *) volk_malloc 
                 (sizeof(float)*
                  (M + cp_len)*3*num_streams,
                  volk_alignment);

  // setting the counters
  sync_index = 0;
  num_samples_processed = 0;
  plateau_start = 0;
  plateau_end = 0;
}

void framesync::print()
{
    printf("ofdmframegen:\n");
    printf("    num subcarriers     :   %-u\n", M);
    printf("      - NULL            :   %-u\n", M_null);
    printf("      - pilot           :   %-u\n", M_pilot);
    printf("      - data            :   %-u\n", M_data);
    printf("    cyclic prefix len   :   %-u\n", cp_len);
    printf("    ");
    ofdmframe_print_sctype(p, M);
}

unsigned long int framesync::get_sync_index()
{
  return sync_index;
}

void framesync::get_G(gr_complex *** _G)
{
  memmove(_G, G, sizeof(gr_complex[M][num_streams][num_streams]));
}

void framesync::reset()
{
  state = STATE_SEEK_PLATEAU;
}

unsigned long long int framesync::get_num_samples_processed()
{
  return num_samples_processed;
}

framesync_states_t
framesync::execute(std::vector<std::complex<float>*>
                   const &in_buff,
                   unsigned int num_samples)
{
  gr_complex x[num_streams];
  if(state == STATE_WAIT)
    state = STATE_MIMO;

  // TODO process a vector of inputs, rather than
  // one by one in a for loop.
  for(unsigned int i = 0; i < num_samples; i++){
    for(unsigned int stream = 0; stream < num_streams; stream++) {
      x[stream] = in_buff[stream][i];
    }

    // FIXME correct frequency offset
    switch(state) {
      case STATE_SEEK_PLATEAU:
	execute_sc_sync(x);
        break;
      case STATE_SAVE_ACCESS_CODES:
	execute_save_access_codes(x);
      case STATE_WAIT:
	break;
      case STATE_MIMO:
	execute_mimo_decode(x);
        break;
      default:
        printf("framesync: state unknown, exiting\n");
        exit(1);
    }
    num_samples_processed++;
  }
  increment_state(); // for testing only
  return state;
}

// for testing only
void framesync::increment_state()
{
  switch(state) {
    case STATE_SEEK_PLATEAU:
      state = STATE_SAVE_ACCESS_CODES;
      break;
    case STATE_SAVE_ACCESS_CODES:
      state = STATE_WAIT;
      break;
    case STATE_WAIT:
      state = STATE_MIMO;
      break;
    case STATE_MIMO:
      state = STATE_MIMO;
      break;
  }
}

void framesync::execute_mimo_decode(gr_complex _x[])
{
}

void framesync::execute_sc_sync(gr_complex _x[])
{
}

void framesync::execute_save_access_codes(gr_complex _x[])
{
}

void framesync::estimate_channel() {
  for(unsigned int i = 0; i < M; i++)
    memmove((((gr_complex (*)[num_streams][num_streams])G) + i),
	    CHANNEL_I, sizeof(gr_complex[num_streams][num_streams]));
}

unsigned long int
framesync::get_plateau_start()
{
  return plateau_start;
}

unsigned long int
framesync::get_plateau_end()
{
  return plateau_end;
}

framesync::~framesync()
{
 for(unsigned int i = 0; i < num_streams; i++)
  {
    free(S1[i]);
    free(s1[i]);
    fftwf_free(X[i]);
    fftwf_free(x[i]);
    fftwf_destroy_plan(fft[i]);
  }
  free(p);
  free(G);
  free(S0);
  free(s0);
  volk_free(volk_buff_fc1);
  volk_free(volk_buff_fc2);
  volk_free(volk_buff_fc3);
  volk_free(volk_buff_f1);
}

void ofdmframe_init_default_sctype(unsigned char * _p, unsigned int _M)
{
    // validate input
    if (_M < 6) {
        fprintf(stderr,"warning: ofdmframe_init_default_sctype(), less than 4 subcarriers\n");
    }

    unsigned int i;
    unsigned int M2 = _M/2;

    // compute guard band
    unsigned int G = 0;
    if (ADD_NULL_CARRIERS) {
      G = _M / 10;
      if (G < 2) G = 2;
    }

    // designate pilot spacing
    unsigned int P = (_M > 34) ? 8 : 4;
    unsigned int P2 = P/2;

    // initialize as NULL
    for (i=0; i<_M; i++)
        _p[i] = OFDMFRAME_SCTYPE_NULL;

    // upper band
    for (i=1; i<M2-G; i++) {
        if ( ((i+P2)%P) == 0 )
            _p[i] = OFDMFRAME_SCTYPE_PILOT;
        else
            _p[i] = OFDMFRAME_SCTYPE_DATA;
    }

    // lower band
    for (i=1; i<M2-G; i++) {
        unsigned int k = _M - i;
        if ( ((i+P2)%P) == 0 )
            _p[k] = OFDMFRAME_SCTYPE_PILOT;
        else
            _p[k] = OFDMFRAME_SCTYPE_DATA;
    }
}

void ofdmframe_validate_sctype(const unsigned char * _p,
                               unsigned int _M,
                               unsigned int * _M_null,
                               unsigned int * _M_pilot,
                               unsigned int * _M_data)
{
    // clear counters
    unsigned int M_null  = 0;
    unsigned int M_pilot = 0;
    unsigned int M_data  = 0;

    unsigned int i;
    for (i=0; i<_M; i++) {
        // update appropriate counters
        if (_p[i] == OFDMFRAME_SCTYPE_NULL)
            M_null++;
        else if (_p[i] == OFDMFRAME_SCTYPE_PILOT)
            M_pilot++;
        else if (_p[i] == OFDMFRAME_SCTYPE_DATA)
            M_data++;
        else {
            fprintf(stderr,"error: ofdmframe_validate_sctype(), invalid subcarrier type (%u)\n", _p[i]);
            exit(1);
        }
    }

    // set outputs
    *_M_null  = M_null;
    *_M_pilot = M_pilot;
    *_M_data  = M_data;
}

void ofdmframe_print_sctype(const unsigned char * _p, unsigned int _M)
{
    unsigned int i;

    printf("[");
    for (i=0; i<_M; i++) {
        unsigned int k = (i + _M/2) % _M;

        switch (_p[k]) {
        case OFDMFRAME_SCTYPE_NULL:     printf(".");    break;
        case OFDMFRAME_SCTYPE_PILOT:    printf("|");    break;
        case OFDMFRAME_SCTYPE_DATA:     printf("+");    break;
        default:
            fprintf(stderr,"error: ofdmframe_print_default_sctype(), invalid subcarrier type\n");
            exit(1);
        }
    }

    printf("]\n");
}

void ofdmframe_init_S0(const unsigned char * _p,
                       unsigned int _M,
                       std::complex<float> * _S0,
                       std::complex<float> * _s0,
                       msequence ms)
{
    unsigned int i;

    unsigned int s;
    unsigned int M_S0 = 0;

    // short sequence
    for (i=0; i<_M; i++) {
        // generate symbol
        s = msequence_generate_symbol(ms,3) & 0x01;

        if (_p[i] == OFDMFRAME_SCTYPE_NULL) {
            // NULL subcarrier
            _S0[i] = 0.0f;
        } else {
            if ( (i%2) == 0 ) {
                // even subcarrer
                _S0[i] = s ? 1.0f : -1.0f;
                M_S0++;
            } else {
                // odd subcarrer (ignore)
                _S0[i] = 0.0f;
            }
        }
    }

    // ensure at least one subcarrier was enabled
    if (M_S0 == 0) {
        fprintf(stderr,"error: ofdmframe_init_S0(), no subcarriers enabled; check allocation\n");
        exit(1);
    }

    // run inverse fft to get time-domain sequence
    // TODO make in independent of liquid
    fft_run(_M, _S0, _s0, LIQUID_FFT_BACKWARD, 0);

    // normalize time-domain sequence level
    // TODO Do this with volk
    float g = 1.0f / sqrtf(M_S0);
    for (i=0; i<_M; i++)
        _s0[i] *= g;
}

void ofdmframe_init_S1(const unsigned char * _p,
                       unsigned int _M,
                       unsigned int _num_access_codes,
                       std::complex<float> * _S1,
                       std::complex<float> * _s1,
                       msequence ms)
{
    unsigned int i, j;

    unsigned int s;
    unsigned int M_S1 = 0;

    // long sequence
    for(j = 0; j < _num_access_codes; j++) {
      for (i=0; i<_M; i++) {
        // generate symbol
        //s = msequence_generate_symbol(ms,1);
        s = msequence_generate_symbol(ms,3) & 0x01;

        if (_p[i] == OFDMFRAME_SCTYPE_NULL) {
            // NULL subcarrier
            (_S1 + _M*j)[i] = 0.0f;
        } else {
            (_S1 + _M*j)[i] = s ? 1.0f : -1.0f;
            M_S1++;
        }
      }
      // run inverse fft to get time-domain sequence
      fft_run(_M, _S1 + _M*j, _s1 + _M*j, LIQUID_FFT_BACKWARD, 0);
    }

    // normalize time-domain sequence level
    float g = 1.0f / sqrtf(M_S1/_num_access_codes);
    for (i=0; i<_M*_num_access_codes; i++)
        _s1[i] *= g;
}
