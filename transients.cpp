#include <cmath>
#include <cstdlib>
#include <iostream>

const float pi = std::acos(-1);
const float root2 = std::sqrt(2);

const int samplerate = 1000;
const int n_bands = 4;
const float min_band = 1.0;
const float max_band = 499.0;
const float transient_threshold = 0.8;

const float fast_decay = std::pow(0.1, 1.0 / samplerate);
const float slow_decay = std::pow(0.5, 1.0 / samplerate);
const float smoothing_decay = std::pow(0.03, 1.0 / samplerate);
const float normalization_decay = std::pow(0.8, 1.0 / samplerate);

// 2nd-order Butterworth lowpass filter
struct BWLowpass {
    float b0, a1, a2, x1, x2, y1, y2;

    BWLowpass() {}

    BWLowpass(float cutoff) {
        // Calculate coefficients
        float k = std::tan(pi * cutoff / samplerate);
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

int main() {
    // Initialize lowpass filters
    BWLowpass filters[n_bands + 1];
    for (int i = 0; i < n_bands + 1; i++) {
        // Calculate where to split the frequency bands
        const float band_exp = std::pow(max_band / min_band, 1.0 / n_bands);
        float cutoff_freq = min_band * std::pow(band_exp, i);
        filters[i] = BWLowpass(cutoff_freq);
    }

    // Initialize peaks
    float fast_peak[n_bands];
    float slow_peak[n_bands];
    for (int i = 0; i < n_bands; i++) {
        fast_peak[i] = 0.0;
        slow_peak[i] = 0.0;
    }

    float smoothing_peak = 0.01;
    float normalization_peak = 0.01;

    // Main loop
    while (true) {
        // TODO: Get sample from microphone
        float sample = (std::rand() % 1000) / 500.0 - 1.0;

        // Calculate DC offset for each frequency band
        float lowpass[n_bands + 1];
        for (int i = 0; i < n_bands + 1; i++) {
            lowpass[i] = filters[i].next(sample);
        }
        float bands[n_bands];
        for (int i = 0; i < n_bands; i++) {
            bands[i] = lowpass[i + 1] - lowpass[i];
        }

        // Calculate transient (formula: sum(fast_peak)^2 - threshold * sum(slow_peak)^2)
        float transient = 0.0;
        for (int i = 0; i < n_bands; i++) {
            fast_peak[i] = std::max(sample, fast_decay * fast_peak[i]);
            slow_peak[i] = std::max(sample, slow_decay * slow_peak[i]);

            transient += fast_peak[i] * fast_peak[i];
            transient -= transient_threshold * slow_peak[i] * slow_peak[i];
        }
        transient = std::max(transient, 0.0f);

        // Apply smoothing and normalization to transients
        smoothing_peak = std::max(transient, smoothing_decay * smoothing_peak);
        normalization_peak = std::max(transient, normalization_decay * normalization_peak);

        float output = smoothing_peak / normalization_peak;
        std::cout << output << std::endl;

        // TODO:
        // 1) Compute how much time was spent during this iteration of the loop
        // 2) Pause execution for (time - 1 / samplerate) seconds
    }
}

