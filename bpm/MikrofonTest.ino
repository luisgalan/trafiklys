// #include <cmath>
// #include <cstdlib>
// #include <iostream>

#include <math.h>
#include <fix_fft.h>

const float pi = acos(-1);
const float root2 = sqrt(2);

const int samplerate = 512;
const long long int sample_interval = 1000000 / samplerate; // microseconds
const int n_bands = 3;
const float min_band = 20.0;
const float max_band = 255.0;
const float noise_floor[n_bands] = { 0.0, 0.0, 0.0 };

const float decay = pow(0.1, 1.0 / samplerate);
const float normalization_decay = pow(0.8, 1.0 / samplerate);
const int buffer_size = 4 * samplerate;

// 2nd-order Butterworth lowpass filter
struct BWLowpass {
    float b0, a1, a2, x1, x2, y1, y2;

    BWLowpass() {}

    BWLowpass(float cutoff) {
        // Calculate coefficients
        float k = tan(pi * cutoff / samplerate);
        float k2 = k * k;
        float s = root2 * k + k2;
        b0 = k2 / (1 + s);
        a1 = 2 * (k2 - 1) / (1 + s);
        a2 = (1 - s) / (1 + s);

        // Initialize state buffer
        x1 = x2 = y1 = y2 = 0.0;
    }

    float next(float sample) {
        // Calculate output sample y0
        float x0 = sample;
        float y0 = b0 * (x0 + 2 * x1 + x2) - a1 * y1 - a2 * y2;

        // Shift state buffer
        x2 = x1;
        x1 = x0;
        y2 = y1;
        y1 = y0;

        return y0;
    }
};


int micPin = A6;
int ledPin = 3;
int val = 0;

BWLowpass filters[n_bands + 1];
float peak[n_bands];

char transient_buffer[buffer_size];
char transient_buffer_im[buffer_size];
int transient_buffer_pos = 0;

void setup() {
  pinMode(micPin, INPUT);
  pinMode(ledPin, OUTPUT);
  Serial.begin(115200);

    // Initialize lowpass filters
    for (int i = 0; i < n_bands + 1; i++) {
        // Calculate where to split the frequency bands
        const float band_exp = pow(max_band / min_band, 1.0 / n_bands);
        float cutoff_freq = min_band * pow(band_exp, i);
        filters[i] = BWLowpass(cutoff_freq);
    }

    // Initialize peaks
    for (int i = 0; i < n_bands; i++) {
        peak[i] = 0.0;
    }
}

// IDE:
// Lav en int8 buffer med længde 4 * samplerate (samplerate = 512)
// Hvis bufferen ikke er fyldt op:
//    - Tilføj transient info til bufferen
// Hvis bufferen er fyldt op:
//    - Vent til det rigtige tidspunkt (når lyset ikke blinker)
//    - Lav FFT over bufferen og estimer bpm

float process_transient(float sample) {
  // Calculate DC offset for each frequency band
  float lowpass[n_bands + 1];
  for (int i = 0; i < n_bands + 1; i++) {
    lowpass[i] = filters[i].next(sample);
  }
  float bands[n_bands];
  for (int i = 0; i < n_bands; i++) {
      bands[i] = lowpass[i + 1] - lowpass[i];
  }

  // Calculate transient
  float transient = 0.0;
  for (int i = 0; i < n_bands; i++) {
      peak[i] = max(abs(sample) - noise_floor[i], decay * peak[i]);
      transient += peak[i] * peak[i];
  }

  normalization_peak = max(transient, normalization_decay * normalization_peak);

  return transient / normalization_peak;
}


// int iteration = 0;

void loop() {
  long long int loop_start = micros();

  float sample = analogRead(micPin);

  float transient = process_transient(sample);
  char fixed_transient = (char)min(max(256 * transient - 128, -128), 127);
  transient_buffer[transient_buffer_pos] = fixed_transient;
  transient_buffer_pos++;

  if (transient_buffer_pos >= buffer_size) {
    transient_buffer_pos = 0;
    // TODO: Perform FFT

    // Initialize imaginary part to zero
    for (int i = 0; i < buffer_size; i++) {
      transient_buffer_im[i] = 0;
    }

    // Perform FFT
    fix_fft(transient_buffer, transient_buffer_im, 11, false);

    // Simple BPM approximation
    int fft_idx = 0;
    float max_magnitude2 = 0.0;

    // Assume bpm is between 80 and 180
    int min_idx = (int)(80.0 / 60 * buffer_size / samplerate)
    int max_idx = (int)(180.0 / 60 * buffer_size / samplerate)

    for (int i = min_idx; i < max_idx; i++) {
      float re = transient_buffer[i];
      float im = transient_buffer_im[i];
      float magnitude2 = re * re + im * im;
      if (magnitude2 > max_magnitude2) {
        fft_idx = i;
        max_magnitude2 = magnitude2;
      }
    }

    float bpm = fft_idx * 60.0 * samplerate / buffer_size;
    Serial.println(bpm);

    // TODO: Interpolate FFT for more accurate bpm
  } else {
    // Wait for next sample interval
    long long int elapsed = micros() - loop_start;
    delayMicroseconds(sample_interval - elapsed);
  }


  // // digitalWrite(ledPin, output>0.7);
  // if (iteration > samplerate / 10) {
  //   iteration = 0;
  //   // Serial.println(output);
  //   // digitalWrite(ledPin, output>0.95);
  //   analogWrite(ledPin, pow(output, 32)*255);
  // }
  // iteration++;

  // if (delay < 0) {
  //   Serial.println("TOO SLOW");
  // }
  // assert(delay >= 0);

  // std::cout << output << std::endl;

  // TODO:
  // 1) Compute how much time was spent during this iteration of the loop
  // 2) Pause execution for (1 / samplerate - time) seconds
}
