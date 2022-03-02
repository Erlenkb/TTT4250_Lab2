import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file_list = ["Free_Field_50cm_Height_d1m.etx", "Free_Field_45cm_Height_d1m.etx", "Imp_tube_12.etx", "Imp_tube_21.etx"]


########## Globl parameters ##############
temperature = 18
T = 273.15 + temperature
fs = 44100
plot_length_seconds = 0.008
p_0 = 101.325 # kPa
rho_0 = 1.186 # kg / m^3
T_0 = 293
x1 = 18

#############################################

def _read_csv(file):
    File_Data = np.loadtxt(file,dtype=float, skiprows=22, max_rows=524287)
    #print(File_Data)
    return File_Data

def _get_c_0(T):
    return 343.2 * np.sqrt(T / 292)

def _getP_0(p_a):
    return rho_0 * ((p_a * T_0) / (p_0 * T))

def _getTransferFunction(p1, p2):
    return p2 / p1

def _getFFT(data1, data2):
    sp = np.pad(data1, (0, 512 - len(data1)), "constant")
    sp2 = np.pad(data2, (0, 512 - len(data2)), "constant")
    sp = np.fft.fft(data1, 512)
    sp = np.trim_zeros(sp, trim="fb")
    sp2 = np.fft.fft(data2, 512)
    sp2 = np.trim_zeros(sp2, trim="fb")
    freq = np.fft.fftfreq(n=len(sp), d=1 / fs)

    return freq,sp, sp2

def _splitDataToLength(array1, array2):
    data1 = array1[1][:int(fs * plot_length_seconds)]
    data1_2 = array1[2][:int(fs * plot_length_seconds)]
    data2 = array2[1][:int(fs * plot_length_seconds)]
    data2_2 = array2[1][:int(fs * plot_length_seconds)]
    return data1, data1_2, data2, data2_2

def _plot_Transfer_Function(freq, H, sp, sp2):

    plt.subplot(1,2,1)
    plt.plot(freq, np.abs(H), label="|H_12|")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude")
    plt.xscale("log")
    plt.grid(which="major")
    plt.grid(which="minor", linestyle=":")
    plt.xscale("log")
    plt.xlim(100, 2000)
    plt.legend()
    plt.ylim(0, 5)
    plt.subplot(1,2,2)
    plt.legend()
    plt.plot(freq, np.abs(sp), label="|p1|")
    plt.plot(freq, np.abs(sp2), label="|p2|")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude")
    plt.xscale("log")
    plt.grid(which="major")
    plt.grid(which="minor", linestyle=":")
    plt.xscale("log")
    plt.xlim(100, 2000)
    # plt.ylim(15, 50)
    plt.legend()

    plt.show()

def _getReflected(freq,c):
    H_R = []
    H_I = []
    for i in freq:
        H_R.append(np.exp(1j*2*np.pi*i / c))
        H_I.append(np.exp(-1j*2*np.pi*i / c))
    return np.array(H_R), np.array(H_I)


def _getReflection(H12, H_I, H_R, freq, c):
    temp = []
    num = H12 - H_I
    den = H_R - H12
    for i in freq:
        temp.append(np.exp((2j*2*np.pi*i*x1) / c))
    temp = np.array(temp)
    num2 = num * temp
    return num2 / den

def _plot_ImpulseResponse_and_FFT(data1, data1_2, data2, data2_2, time, f_upper):

    plt.subplot(1, 2, 1)
    plt.plot(time, data1, color="blue", label="Microphone A")
    plt.plot(time, data1_2, color="red", label="Microphone B")
    plt.grid()
    plt.ylabel("Amplitude [Pa]")
    plt.xlabel("Time [ms]")
    plt.legend()

    plt.subplot(1, 2, 2)
    freq, sp, sp2 = _getFFT(data1, data1_2)

    plt.plot(np.fft.fftshift(freq), np.fft.fftshift(np.abs(sp)), color="blue", label="Microphone A")
    plt.plot(np.fft.fftshift(freq), np.fft.fftshift(np.abs(sp2)), color="red", label="Microphone B")

    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude")
    plt.xscale("log")
    plt.grid(which="major")
    plt.grid(which="minor", linestyle=":")
    plt.xscale("log")
    #plt.xlim(100, 2000)
    #plt.ylim(15, 50)
    plt.legend()

    plt.show()

Imp2_1_data = np.transpose(_read_csv(file_list[3]))
Imp_1_2_data = np.transpose(_read_csv(file_list[2]))
time_length = Imp_1_2_data.shape[1]
time_axis = np.linspace(0, 8, int(fs * plot_length_seconds))

c = _get_c_0(T)

data1, data1_2, data2, data2_2 = _splitDataToLength(Imp_1_2_data, Imp2_1_data)

freq, sp, sp2 = _getFFT(data1,data1_2)

freq = np.fft.fftshift(freq)
sp = np.fft.fftshift(sp)
sp2 = np.fft.fftshift(sp2)

H12 = _getTransferFunction(sp, sp2)

_plot_Transfer_Function(freq,H12,sp,sp2)

_plot_ImpulseResponse_and_FFT(data1, data1_2, data2, data2_2,time_axis,2000)

H_R, H_I = _getReflected(freq,c)

Reflection_R = _getReflection(H12, H_I, H_R, freq, c)

plt.plot(freq, np.abs(Reflection_R), label="Reflection Coefficient R")

plt.xlabel("Frequency [Hz]")
plt.ylabel("Magnitude")
plt.xscale("log")
plt.grid(which="major")
plt.grid(which="minor", linestyle=":")
plt.xscale("log")
plt.xlim(100, 2000)

plt.legend()

plt.show()




