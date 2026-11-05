# One FMCW frame of complex ADC samples (dechirped IQ).
std_msgs/Header header
uint32 num_rx
uint32 num_chirps
uint32 num_samples
# Interleaved real/imaginary float32, layout [rx][chirp][sample].
float32[] iq
