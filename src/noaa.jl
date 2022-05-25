using WAV
using DSP
using FFTW
using Images

using Plots
backend(:plotly)
gr(size = (800, 300), legend = true,  background_color = RGB(0.2, 0.2, 0.2),)

file = "/home/konrad/Downloads/APT-testfiles/176400_clear.wav"

frame_len = 4160


# read the data and samplerate from the file 
data, sampleRate = wavread(file)


data = DSP.Filters.resample(data, float(frame_len * 2) / float(sampleRate))
# demodulate the signal 
data = DSP.Util.hilbert(data[1:1:end])


lines = 320
start = 400


length_sum = ((4160 * 2) * lines)
formated = abs.(output_data[start:1:length_sum + start - 1]) 
print(length(formated))

formated = abs.(reshape(formated, 4160 * 2, lines))


gui(heatmap(formated, c=:grayscale, transpose = true))