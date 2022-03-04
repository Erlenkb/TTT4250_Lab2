import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

file_list = ["Free_Field_50cm_Height_d1m.etx", "Free_Field_45cm_Height_d1m.etx", "Imp_tube_12.etx", "Imp_tube_21.etx"]


########## Globl parameters ##############
temperature = 18
T = 273.15 + temperature
fs = 44100
plot_length_seconds = 0.1
p_0 = 101.325 # kPa
rho_0 = 1.186 # kg / m^3
T_0 = 293
x1 = 0.18
s = 0.08

#############################################

def _nextpow2(i):
    n = 1
    while n < i: n *= 2
    return n



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
    sp = np.pad(data1, (0, _nextpow2(len(data1)) - len(data1)), "constant")
    sp2 = np.pad(data2, (0, _nextpow2(len(data2)) - len(data2)), "constant")
    sp = np.fft.fft(sp, _nextpow2(len(data1)))
    sp = np.trim_zeros(sp, trim="fb")
    sp2 = np.fft.fft(sp2, _nextpow2(len(data2)))
    sp2 = np.trim_zeros(sp2, trim="fb")
    freq = np.fft.fftfreq(n=len(sp), d=1 / fs)

    return freq, sp, sp2

def _splitDataToLength(array1, array2):
    data1 = array1[1][:int(fs * plot_length_seconds)]
    data1_2 = array1[2][:int(fs * plot_length_seconds)]
    data2 = array2[1][:int(fs * plot_length_seconds)]
    data2_2 = array2[1][:int(fs * plot_length_seconds)]
    return data1, data1_2, data2, data2_2

def _plot_Transfer_Function(freq, H, sp, sp2):
    fig = plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.semilogx(freq, np.abs(H), label="|H_12|")
    ax1.grid(which="major")
    ax1.grid(which="minor", linestyle=":")
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_ylabel("Amplitude [dB]")
    ax1.set_title(" Burde kanskje finne en tittel")

    ax2.semilogx(freq, 20*np.log10(sp), label="|p1|")
    ax2.semilogx(freq, 20 * np.log10(sp2), label="|p2|")
    ax2.grid(which="major")
    ax2.grid(which="minor", linestyle=":")
    ax2.set_xlabel("Frequency [Hz]")
    ax2.set_ylabel("Amplitude [dB]")
    ax2.set_title(" Burde kanskje finne en tittel")

    ax1.set_xlim(100,2000)
    ax1.set_ylim(0,3)
    ax2.set_xlim(100,2000)
    ax2.set_ylim(10,50)

    ax1.legend()
    ax2.legend()
    plt.show()
    fig.savefig("Transfer_Function&FFT.png")
    plt.close(fig)

def _getReflected(freq,c):
    H_R = []
    H_I = []

    Angular_freq = 2*np.pi*freq

    for i in Angular_freq:
        H_R.append(np.exp(s*1j*i / c))
        H_I.append(np.exp(s*-1j*i / c))
    return np.array(H_R), np.array(H_I)


def _getReflection(H12, H_I, H_R, freq, c):
    temp = []
    num = H12 - H_I
    den = H_R - H12
    Angular_freq = 2*np.pi*freq
    for i in Angular_freq:
        temp.append(np.exp((2j*2*np.pi*i*x1) / c))
    temp = np.array(temp)
    num2 = num * temp
    return num2 / den


def _plot_ImpulseResponse_and_FFT(data1, data1_2, data2, data2_2, time, f_upper):
    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)

    ax1.plot(time, data1, color="blue", label="Microphone A")
    ax1.plot(time, data1_2, color="red", label="Microphone B")
    ax1.grid(which="major")
    #ax1.grid(which="minor", linestyle=":")
    ax1.set_ylabel("Amplitude [Pa]")
    ax1.set_xlabel("Time [ms]")
    ax1.legend()


    freq, sp, sp2 = _getFFT(data1, data1_2)

    ax2.semilogx(np.fft.fftshift(freq), 20*np.log10(np.fft.fftshift(sp)), color="blue", label="Microphone A")
    ax2.semilogx(np.fft.fftshift(freq), 20*np.log10(np.fft.fftshift(sp2)), color="red", label="Microphone B")
    ax2.grid(which="major")
    ax2.grid(which="minor", linestyle=":")
    ax2.set_ylabel("Magnitude [|Pa|]")
    ax2.set_xlabel("Frequency [Hz]")
    ax2.legend()
    ax2.set_xlim(100, 2000)
    ax2.set_ylim(10,50)

    ax2.legend()
    fig.savefig("Impulse_Response&FFT_First_Microphone_setup.png")
    plt.show()
    plt.close(fig)



def _plot_Reflection_Coefficient(R,freq):
    fig = plt.figure(figsize=(7,4))
    plt.plot(freq, np.abs(R), label="Reflection Coefficient R")
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Magnitude")
    plt.xscale("log")
    plt.grid(which="major")
    plt.grid(which="minor", linestyle=":")
    plt.xscale("log")
    plt.xlim(100, 2000)
    plt.ylim(0,1)
    plt.legend()
    plt.show()
    fig.savefig("Reflection_Coefficient.png")
    plt.close(fig)

def _plot_AbsorptionCoefficient_and_Impedance(R, freq):
    alpha = 1 - (np.abs(R))**2
    Z = (1+R) / (1-R)

    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.semilogx(freq,alpha, label="Absorption coefficient")
    ax1.grid(which="major")
    ax1.grid(which="minor", linestyle=":")
    ax1.set_ylabel("Value")
    ax1.set_xlabel("Frequency [Hz]")
    ax1.set_xlim(100,2000)
    ax1.set_ylim(0,1)
    ax1.legend()

    #ax2.semilogx(freq, Z.real, label="Re[Z]")
    #ax2.semilogx(freq, Z.imag, label="Im[Z]")
    ax2.semilogx(freq, np.abs(Z), label="|Z|")
    ax2.grid(which="major")
    ax2.grid(which="minor", linestyle=":")
    ax2.set_ylabel("Value")
    ax2.set_xlabel("Frequency [Hz]")
    ax2.legend()
    ax2.set_xlim(100,2000)
    ax2.set_ylim(-10,10)
    plt.show()
    fig.savefig("Alpha&Impedance.png")
    plt.close(fig)

def _calculate_Corrected_TrensferFunc(H1, H2):
    HC = np.sqrt(H1 / H2)
    return H1 / HC


#### The boring calculations####################

Imp2_1_data = np.transpose(_read_csv(file_list[3]))
Imp_1_2_data = np.transpose(_read_csv(file_list[2]))
time_length = Imp_1_2_data.shape[1]
time_axis = np.linspace(0, 8, int(fs * plot_length_seconds))
c = _get_c_0(T)

################################################


###### Fetching the correct data lengths and calculcate the FFT

data1, data1_2, data2, data2_2 = _splitDataToLength(Imp_1_2_data, Imp2_1_data)

freq, sp, sp2 = _getFFT(data1,data1_2)

freq2, sp3, sp4 = _getFFT(data2, data2_2)

#### Shift the FFT and frecvector to prevent the annoying straight line in the plot due to frecvec start at positive then negative.
freq = np.fft.fftshift(freq)
sp = np.fft.fftshift(sp)
sp2 = np.fft.fftshift(sp2)
freq2 = np.fft.fftshift(freq2)
sp3 = np.fft.fftshift(sp3)
sp4 = np.fft.fftshift(sp4)

#############################
H12 = _getTransferFunction(sp,sp2)        #### Should I put the 10*np.log10(FFTvalues) in here? Or just the FFTvalues?
H12_2 = _getTransferFunction(sp4, sp3)
############################

#H12 = _calculate_Corrected_TrensferFunc(H12, H12_2)

_plot_Transfer_Function(freq,H12,sp,sp2)


_plot_ImpulseResponse_and_FFT(data1, data1_2, data2, data2_2,time_axis,2000)

H_R, H_I = _getReflected(freq,c)


#### Create R and plot R, alpha and Z

Reflection_R = _getReflection(H12, H_I, H_R, freq, c)
R = Reflection_R
_plot_Reflection_Coefficient(R, freq)
_plot_AbsorptionCoefficient_and_Impedance(R,freq)







