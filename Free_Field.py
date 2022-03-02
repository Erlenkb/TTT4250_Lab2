import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

file_list = ["Free_Field_50cm_Height_d1m.etx", "Free_Field_45cm_Height_d1m.etx", "Imp_tube_12.etx", "Imp_tube_21.etx"]

file1 = "Free_Field_50cm_Height_d1m.etx"
file2 = "Free_Field_45cm_Height_d1m.etx"

####### Important values #######

temperature = 18
T = 273.15 + temperature
fs = 44100
plot_length_seconds = 0.008

################################



def _read_csv(file):
    File_Data = np.loadtxt(file,dtype=float, skiprows=22, max_rows=524287)
    return File_Data

def _get_max_f_upper(c_0):
    f_upper_s_argument = 0.45 * c_0 / 0.08
    f_upper_d_argument = 0.58 * c_0 / 0.1
    print("f_upper as function of s: {0:.2f}Hz\nf_upper as function of d: {1:.2f}".format(f_upper_s_argument,f_upper_d_argument))
    return 0.45 * c_0 / 8

def _plot_ImpulseResponse_and_FFT(array1, array2, time, f_upper):
    data1 = array1[1][:int(fs * plot_length_seconds)]
    data2 = array2[1][:int(fs * plot_length_seconds)]

    plt.subplot(1, 2, 1)
    plt.plot(time, data1, color="blue", label="Upper microphone")
    plt.plot(time, data2, color="red", label="Lower microphone")
    plt.grid()
    plt.ylabel("Amplitude [Pa]")
    plt.xlabel("Time [ms]")
    plt.legend()

    plt.subplot(1, 2, 2)
    sp = np.pad(data1, (0, 512 - len(data1)), "constant")
    sp2 = np.pad(data2, (0, 512 - len(data2)), "constant")
    sp = np.fft.fft(data1, 512)
    sp = np.trim_zeros(sp, trim="fb")
    sp2 = np.fft.fft(data2, 512)
    sp2 = np.trim_zeros(sp2, trim="fb")

    freq = np.fft.fftfreq(n=len(sp), d=1 / fs)

    plt.plot(np.fft.fftshift(freq), np.fft.fftshift(np.abs(sp)), color="blue", label="Upper microphone")
    plt.plot(np.fft.fftshift(freq), np.fft.fftshift(np.abs(sp2)), color="red", label="Upper microphone")

    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude")
    plt.xscale("log")
    plt.grid(which="major")
    plt.grid(which="minor", linestyle=":")
    plt.xscale("log")
    plt.xlim(82, 2000)
    plt.ylim(15, 50)
    plt.legend()

    plt.show()


if __name__ == '__main__':

    file1_Data = np.transpose(_read_csv(file1))
    file2_Data = np.transpose(_read_csv(file2))
    print(file1_Data.shape)

    c_0 = _get_c_0(T)
    f_upper_s = _get_max_f_upper(c_0)


    time_length = file1_Data.shape[1]
    time_axis = np.linspace(0, 8, int(fs * plot_length_seconds))

    _plot_ImpulseResponse_and_FFT(file1_Data, file2_Data, time_axis, f_upper_s)









