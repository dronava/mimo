/*
 * Copyright (c) 2014 Manu T S
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

#include <math.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <getopt.h>
#include <liquid/liquid.h>
#include <liquid/mimo.h>
#include <ctime>

#include <uhd/utils/thread_priority.hpp>
#include <uhd/usrp/multi_usrp.hpp>
#include <uhd/utils/safe_main.hpp>
#include <gnuradio/filter/fft_filter.h>
#include <gnuradio/gr_complex.h>
#include <volk/volk.h>

#include <boost/thread.hpp>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

unsigned int num_frames_detected;
static bool verbose;

unsigned int fill_tx_sig(std::complex<float> ** tx_sig,
                         unsigned int tx_sig_len,
                         unsigned int k,
                         unsigned int m,
                         float beta,
                         float tone_amp,
                         unsigned int * num_frames,
                         bool save);

int callback(unsigned char *  _header,
             int              _header_valid,
             unsigned char *  _payload,
             unsigned int     _payload_len,
             int              _payload_valid,
             framesyncstats_s _stats,
             void *           _userdata);

int UHD_SAFE_MAIN(int argc, char **argv)
{
  uhd::set_thread_priority_safe();

  // operating parameters
  double cent_freq;         // center frequency of transmission
  double samp_rate;         // usrp samping rate
  float tone_amp;           // tone amplitude
  double txgain;            // tx frontend gain
  double rxgain;            // rx frontend gain
  double num_secs;          // number of seconds to operate

  // default values
  double d_cent_freq = 2600.0e6;
  double d_samp_rate = 1000e3;
  float d_tone_amp   = 0.25;
  double d_txgain    = 20.0;
  double d_rxgain    = 20.0;
  double d_num_secs  = 5.0;
  bool d_verbose     = true;

  //set the operating parameters
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "help message")
    ("freq", po::value<double>(&cent_freq)->default_value(d_cent_freq), "RF center frequency in Hz")
    ("rate", po::value<double>(&samp_rate)->default_value(d_samp_rate), "USRP Sampling rate")
    ("rrc_amplitude", po::value<float>(&tone_amp)->default_value(d_tone_amp), "rrc pulse shape amplitude")
    ("tx_gain", po::value<double>(&txgain)->default_value(d_txgain), "TX Front end gain")
    ("rx_gain", po::value<double>(&rxgain)->default_value(d_rxgain), "RX Front end gain")
    ("duration", po::value<double>(&num_secs)->default_value(d_num_secs), "Number of seconds to run")
    ("verbose", po::value<bool>(&verbose)->default_value(d_verbose), "Verbose")
    ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  //print the help message
  if (vm.count("help")) {
    std::cout << boost::format("mimo-test %s") % desc << std::endl;
    std::cout
      << std::endl
      << "Testing MIMO.\n"
      << std::endl;
      return 0;
  }

  // streamer properties
  std::string cpu = "fc32";           // cpu format for the streamer
  std::string wire = "sc16";          // wire formate for the streamer

  // pulse shape properties
  unsigned int k  = 16;               // samples/symbol
  unsigned int m  = 11;               // filter delay (symbols)
  float beta      = 0.35f;            // excess bandwidth factor

  // channel properties
  size_t num_rx_chans;
  size_t num_tx_chans;
  std::vector<size_t> tx_chans;
  std::vector<size_t> rx_chans;

  // rest
  bool save = true;
  size_t reps;
  size_t tx_buff_len;
  size_t rx_buff_len;
  size_t tx_sig_len;
  size_t rx_sig_len;
  size_t num_rx_samps, num_tx_samps;
  size_t num_samps_sent, num_samps_rcvd;
  float timeout = 0.1;
  unsigned int num_tx_frames;
  void * user_data;

  std::complex<float> * tx_buff[2];
  std::complex<float> * tx_sig[2];
  std::complex<float> * rx_buff[2];
  std::complex<float> * rx_sig[2];

  // setting up the link
  uhd::device_addr_t tx_addr, rx_addr;
  tx_addr["addr0"] = "134.147.118.216";
  rx_addr["addr0"] = "134.147.118.217";
  uhd::usrp::multi_usrp::sptr tx = uhd::usrp::multi_usrp::make(tx_addr);
  uhd::usrp::multi_usrp::sptr rx = uhd::usrp::multi_usrp::make(rx_addr);
  uhd::usrp::subdev_spec_t tx_subdev_spec("A:0 B:0");
  uhd::usrp::subdev_spec_t rx_subdev_spec("A:0 B:0");
  tx->set_tx_subdev_spec(tx_subdev_spec, uhd::usrp::multi_usrp::ALL_MBOARDS);
  rx->set_rx_subdev_spec(rx_subdev_spec, uhd::usrp::multi_usrp::ALL_MBOARDS);

  num_tx_chans = tx->get_tx_num_channels();
  // set freq, rate, gain, antenna.
  std::cout << "Setting rate, gain and freq for TX\n";
  for (size_t chan = 0; chan < num_tx_chans; chan++) {
    std::cout << "Parameters for Channel " << chan <<"\n";
    tx->set_tx_rate(samp_rate, chan);
    std::cout << "TX Rate :";
    std::cout << tx->get_tx_rate(chan) << "\n";
    uhd::tune_request_t tx_tune_request(cent_freq);
    uhd::tune_result_t tx_tune_result;
    tx_tune_result = tx->set_tx_freq(tx_tune_request, chan);
    std::cout << "Transmit Tune Result" << "\n";
    std::cout << tx_tune_result.to_pp_string();
    tx->set_tx_gain(txgain, chan);
    std::cout << "Transmit Gain :" << tx->get_tx_gain(chan) << "dB\n";
    tx->set_tx_antenna("TX/RX", chan);
    std::cout << "Transmit Antenna :" << tx->get_tx_antenna(chan) << "\n";
    std::cout << "Transmit subdevice :" << tx->get_tx_subdev_name(chan) << "\n";
    tx_chans.push_back(chan);
  }

  num_rx_chans = rx->get_rx_num_channels();
  // set freq, rate, gain, antenna.
  std::cout << "Setting rate, gain and freq for RX\n";
  for (size_t chan = 0; chan < num_rx_chans; chan++) {
    std::cout << "Parameters for Channel " << chan <<"\n";
    rx->set_rx_rate(samp_rate, chan);
    std::cout << "RX Rate :";
    std::cout << rx->get_rx_rate(chan) << "\n";
    uhd::tune_request_t rx_tune_request(cent_freq);
    uhd::tune_result_t rx_tune_result;
    rx_tune_result = rx->set_rx_freq(rx_tune_request, chan);
    std::cout << "Receive Tune Result" << "\n";
    std::cout << rx_tune_result.to_pp_string();
    rx->set_rx_gain(rxgain, chan);
    std::cout << "Receive Gain :" << rx->get_rx_gain(chan) << "dB\n";
    rx->set_rx_antenna("TX/RX", chan);
    std::cout << "Receive Antenna :" << rx->get_rx_antenna(chan) << "\n";
    std::cout << "Receive subdevice :" << rx->get_rx_subdev_name(chan) << "\n";
    rx_chans.push_back(chan);
  }

  // create a tx streamer
  uhd::stream_args_t tx_stream_args(cpu, wire);
  tx_stream_args.channels = tx_chans;
  uhd::tx_streamer::sptr tx_stream = tx->get_tx_stream(tx_stream_args);

  // create an rx streamer
  uhd::stream_args_t rx_stream_args(cpu, wire);
  rx_stream_args.channels = rx_chans;
  uhd::rx_streamer::sptr rx_stream = rx->get_rx_stream(rx_stream_args);

  // allocate memory to store the transmit signal
  tx_buff_len = tx_stream->get_max_num_samps();
  tx_sig_len = (size_t)(num_secs*samp_rate);
  std::cout << "Transmit Buffer Length :" << tx_buff_len << "\n";
  for (size_t chan = 0; chan < num_tx_chans; chan++) {
    tx_buff[chan] = (std::complex<float> *)malloc(tx_buff_len*sizeof(std::complex<float>));
    tx_sig[chan] = (std::complex<float> *)malloc(tx_sig_len*sizeof(std::complex<float>));
  }
  // allocate memory to store the received samples
  rx_buff_len = rx_stream->get_max_num_samps();
  rx_sig_len = (size_t)(num_secs*samp_rate);
  std::cout << "Receive Buffer Length :" << rx_buff_len << "\n";
  for (size_t chan = 0; chan < num_rx_chans; chan++) {
    rx_buff[chan] = (std::complex<float> *)malloc(rx_buff_len*sizeof(std::complex<float>));
    rx_sig[chan] = (std::complex<float> *)malloc(rx_sig_len*sizeof(std::complex<float>));
  }

  reps = 1;
  for(size_t run = 0; run < reps; run++) {
    tx_sig_len = fill_tx_sig(tx_sig, tx_sig_len, k, m, beta, tone_amp,
        &num_tx_frames, save);

    tx->set_time_now(uhd::time_spec_t(0.0), uhd::usrp::multi_usrp::ALL_MBOARDS);
    rx->set_time_now(uhd::time_spec_t(0.0), uhd::usrp::multi_usrp::ALL_MBOARDS);

    uhd::stream_cmd_t stream_cmd(uhd::stream_cmd_t::STREAM_MODE_START_CONTINUOUS);
    stream_cmd.stream_now = false;
    stream_cmd.time_spec = uhd::time_spec_t(1.0);
    rx_stream->issue_stream_cmd(stream_cmd);

    uhd::tx_metadata_t txmd;
    txmd.start_of_burst = true;
    txmd.end_of_burst = false;
    uhd::rx_metadata_t rxmd;
    num_tx_samps = 0;
    num_rx_samps = 0;

    timeout = 1.1;
    while(true)
    {
      // populate the transmit buffer
      for (size_t chan = 0; chan < num_tx_chans; chan++) {
        memmove(tx_buff[chan],
                tx_sig[chan] + num_tx_samps,
                sizeof(std::complex<float>)*tx_buff_len);
      }
      // transmit the transmit buffer
      num_samps_sent = tx_stream->send(tx_buff,
                                       tx_buff_len,
                                       txmd);
      // receive samples to the receive-buffer.
      num_samps_rcvd = rx_stream->recv(rx_buff,
                                       rx_buff_len,
                                       rxmd,
                                       timeout);
      if (num_samps_sent != tx_buff_len) {
        std::cout << "\nRequested " << tx_buff_len << " samples to be sent\n";
        std::cout << "send returned with " << num_samps_sent << " samples\n";
      }
      if (num_samps_rcvd != rx_buff_len) {
        std::cout << "\nRequested " << rx_buff_len << " samples to be read\n";
        std::cout << "recv returned with " << num_samps_rcvd << " samples\n";
      }
      if(rxmd.error_code) {
        std::cerr << "Receive Stream Error Code :" << rxmd.error_code << "\n";
        break;
      }
      if(num_rx_samps + num_samps_rcvd > rx_sig_len)
        break;
      if(num_tx_samps + tx_buff_len > tx_sig_len)
        break;
      // save the received samples
      for (size_t chan = 0; chan < num_rx_chans; chan++)
        memmove(rx_sig[chan] + num_rx_samps,
                rx_buff[chan],
                sizeof(std::complex<float>)*num_samps_rcvd);
      num_rx_samps += num_samps_rcvd;
      num_tx_samps += num_samps_sent;
      txmd.start_of_burst = false;
      timeout = 0.1;
    }
    // send the last packet
    txmd.end_of_burst = true;
    tx_stream->send(tx_buff,
                    0,
                    txmd);
    // stop rx streaming
    stream_cmd.stream_mode = uhd::stream_cmd_t::STREAM_MODE_STOP_CONTINUOUS;
    rx_stream->issue_stream_cmd(stream_cmd);
  }
  if(save){
    FILE * f_rx_sig1;
    FILE * f_rx_sig2;
    f_rx_sig1 = fopen("/tmp/rx_sig1", "wb");
    f_rx_sig2 = fopen("/tmp/rx_sig2", "wb");
    fwrite((void *)(rx_sig[0]), sizeof(std::complex<float>), num_rx_samps, f_rx_sig1);
    fwrite((void *)(rx_sig[1]), sizeof(std::complex<float>), num_rx_samps, f_rx_sig2);
    fclose(f_rx_sig1);
    fclose(f_rx_sig2);
  }

  user_data = &samp_rate;
  liquid::mimo::framesync fs(callback, user_data, k, m, beta);
  fs.work(rx_sig[0], num_rx_samps);
  std::cout << "Number of PN1 received = "
    << fs.get_frame1_count() << std::endl;
  std::cout << "Number of PN2 received = "
    << fs.get_frame2_count() << std::endl;
  std::cout << "Number of PN3 received = "
    << fs.get_frame3_count() << std::endl;
  std::cout << "tx_sig_len = "
    << tx_sig_len << std::endl;
  std::cout << "num_tx_frames = "
    << num_tx_frames << std::endl;

  for (size_t chan = 0; chan < num_tx_chans; chan++) {
    free(tx_sig[chan]);
    free(tx_buff[chan]);
  }
  for (size_t chan = 0; chan < num_rx_chans; chan++) {
    free(rx_sig[chan]);
    free(rx_buff[chan]);
  }
  std::cout << "num_tx_samps = " << num_tx_samps << std::endl;
  std::cout << "num_rx_samps = " << num_rx_samps << std::endl;
  return EXIT_SUCCESS;
}

unsigned int fill_tx_sig(std::complex<float> ** tx_sig,
                         unsigned int tx_sig_len,
                         unsigned int k,
                         unsigned int m,
                         float beta,
                         float tone_amp,
                         unsigned int * num_frames,
                         bool save) {
  liquid::mimo::framegen packer(k, m, beta);
  packer.set_gains(tone_amp, tone_amp);
  unsigned int num_samps = packer.work(tx_sig, tx_sig_len);
  *num_frames = packer.get_frame_count();
  if(save){
    FILE * f_tx_sig1;
    FILE * f_tx_sig2;
    f_tx_sig1 = fopen("/tmp/tx_sig1", "wb");
    f_tx_sig2 = fopen("/tmp/tx_sig2", "wb");
    fwrite((void *)(tx_sig[0]), sizeof(std::complex<float>), num_samps, f_tx_sig1);
    fwrite((void *)(tx_sig[1]), sizeof(std::complex<float>), num_samps, f_tx_sig2);
    fclose(f_tx_sig1);
    fclose(f_tx_sig2);
  }
  return num_samps;
}

void fill_tx_sig(std::vector<std::complex<float> *> tx_sig,
                 size_t tx_sig_len,
                 std::complex<float> * pn1,
                 std::complex<float> * pn2,
                 size_t seq_len,
                 float tone_amp,
                 firinterp_crcf interp,
                 bool save
                ){
  const std::complex<float> I(0.0, 1.0);
  const std::complex<float> Z(0.0, 0.0);
  size_t tx_sig_counter;
  size_t k = 16;
  bool cont;
  size_t samps_counter;
  std::complex<float> * smb;
  size_t smb_counter;
  smb = (std::complex<float> *)malloc(sizeof(std::complex<float>)*seq_len*2);
  std::complex<float> * samps;
  samps = (std::complex<float> *)malloc(sizeof(std::complex<float>)*k);
  for(unsigned int i = 0; i < seq_len; i++)
  {
    smb[i] = pn1[i]*tone_amp;
    smb[i + seq_len] = I*pn2[i]*tone_amp;
  }

  smb_counter = 0;
  tx_sig_counter = 0;
  cont = true;
  while(cont) {
    firinterp_crcf_execute(interp, smb[smb_counter], samps);
    for(samps_counter = 0; samps_counter < k; samps_counter++)
    {
      tx_sig[0][tx_sig_counter] = samps[samps_counter].real() + Z;
      tx_sig[1][tx_sig_counter] = Z + I*samps[samps_counter].imag();
      tx_sig_counter++;
      if(!(tx_sig_counter < tx_sig_len)) {
        cont = false;
        break;
      }
    }
    smb_counter++;
    smb_counter = smb_counter%(seq_len*2);
  }
  if(save){
    FILE * f_tx_sig1;
    FILE * f_tx_sig2;
    f_tx_sig1 = fopen("/tmp/tx_sig1", "wb");
    f_tx_sig2 = fopen("/tmp/tx_sig2", "wb");
    fwrite((void *)(tx_sig[0]), sizeof(std::complex<float>), tx_sig_len, f_tx_sig1);
    fwrite((void *)(tx_sig[1]), sizeof(std::complex<float>), tx_sig_len, f_tx_sig2);
    fclose(f_tx_sig1);
    fclose(f_tx_sig2);
  }
}

// callback function
int callback(unsigned char *  _header,
             int              _header_valid,
             unsigned char *  _payload,
             unsigned int     _payload_len,
             int              _payload_valid,
             framesyncstats_s _stats,
             void *           _userdata)
{
    if (verbose) {
        // compute true carrier offset
        double samplerate = *((double*)_userdata);
        float cfo = _stats.cfo / M_PI * samplerate * 1e-3f;
        printf("***** rssi=%7.2fdB evm=%7.2fdB, cfo=%7.2fkHz, ", _stats.rssi, _stats.evm, cfo);

        if (_header_valid) {
            unsigned int packet_id = (_header[0] << 8 | _header[1]);
            printf("rx packet id: %6u", packet_id);
            if (_payload_valid) printf("\n");
            else                printf(" PAYLOAD INVALID\n");
        } else {
            printf("HEADER INVALID\n");
        }
    } else {
    }

    // update global counters
    num_frames_detected++;

    return 0;
}
