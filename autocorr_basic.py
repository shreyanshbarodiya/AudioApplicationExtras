from math import *
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
from scipy.io.wavfile import read

# def envelope(data_in):
#     data_out = [data_in[0]]
#     data_index = [0]
#     length = len(data_in)
#     for k in range(1, length):
#         if(data_in[k]>data_in[k-1]):
#             data_out.append(data_in[k])
#             data_index.append(k)
#
#     return data_index, data_out
#
# def envelope_interpolate(data, data_index, data_env):
#     data_out=[data[0]]
#     pos = 0
#     length = len(data)
#     data_tmp = data.tolist()    # converying numpy array to list
#
#     for m in range(1, int(length/10)):
#         k = 10*m
#         if k in data_index == True:
#             data_out.append(data_tmp[k])
#             pos = data_index.index(k)       # position of data in envelope data
#         else:
#             dist = data_index[pos+1]-data_index[pos]
#             wt1 = k - data_index[pos]
#             wt2 = data_index[pos+1] - k
#             data_val = (wt1*data_env[pos+1]+wt2*data_env[pos])/(wt1+wt2)
#             data_out.append(data_val)
#         print k
#
#     return data_out

def envelope_detect(data_in, samp_rate):
    data_tmp = signal.decimate(data_in, 10)     # downsampling
    data_env = abs(signal.hilbert(data_tmp))
    Fc = 50.0/(samp_rate*10)     # cutoff freq of 50Hz
    a = 1
    b = signal.firwin(121, cutoff = Fc, window='hamming')
    print b
    #input('Press any key')

    data_out = signal.lfilter(b,a, data_env)
    return data_out, int(samp_rate/10)

samp_rate, data_in = read('../Data/NST_Normal.wav')
#samp_rate, data_in = read('../Data/fhr-02_mono.wav')
#samp_rate, data_in = read('../Data/Tachycardia.wav')

tmin = 0
data_in = data_in[tmin*samp_rate:(tmin+10)*samp_rate]
time_len = len(data_in)
time_data_in = np.linspace(tmin, time_len/samp_rate+tmin, time_len)
#plt.subplot(311)
plt.plot(time_data_in, data_in, 'b')

data, samp_rate = envelope_detect(data_in, samp_rate)   # we have decimated the signal by a factor of 10
time_len = len(data)
time_data = np.linspace(tmin, time_len/samp_rate+tmin, time_len)
#plt.subplot(312)
plt.plot(time_data, data, 'r')

# data_tmp = data_in[0:500000]
# data_index, data_env = envelope(data_tmp)


# data = envelope_interpolate(data_tmp, data_index, data_env)
#
# plt.plot(data)
plt.show()

win_len = int(1.5*samp_rate)     # Auto-correlation measured over a window of 1.5s
tot_len = len(data)
win_shift = 0.2*samp_rate  # Window shifted by 200ms for nex ACF calculation

num_iter = int((tot_len - win_len)/win_shift) + 1
#rate = []
#time = []
rate = np.zeros(num_iter)
time = np.linspace(0, (num_iter-1)*win_shift, num_iter)/samp_rate + tmin
last = 0

k = 0
# traversing the entire data

for i in range(0, num_iter):
    win_bot = i*win_shift
    win_data = data[win_bot:win_bot+win_len]

    acf = []
    fac = 10
    min_period = 60.0/300  # Max bpm = 300
    max_period = 60.0/60   # Min bpm = 60
    acf_min = int(samp_rate*min_period/fac)
    acf_max = int(samp_rate*max_period/fac)


    # determining acf values over window
    for j in range(acf_min, acf_max):          # ACF is max at zero, so need to exclude that
        win_data_1 = win_data[0:win_len-j*fac]
        win_data_2 = win_data[j*fac:win_len]
        acf_val = sum(win_data_1*win_data_2)
        acf.append(acf_val)

    print len(range(acf_min, acf_max)), len(acf)
# determining period from acf
    period = (acf_min + acf.index(max(acf)))*fac    # period in terms of samples
    print "Period is: ", period, " at time ", i*0.2
    plt.plot(range(acf_min, acf_max), acf)
    plt.show()
    plt.title('ACF')
    per_time = float(period*1000)/samp_rate     # period in ms
    rate_tmp = 60000/per_time
    print "Calculated rate: ", rate_tmp, " at time ", i*0.2

    if i==0 or abs(rate_tmp - rate[last]) <= 50:
        #rate.append(rate_tmp)
        #time.append(i*win_shift/samp_rate)
        rate[i] = rate_tmp
        last = i
    else:
        rate[i] = np.nan

    k = k+1
    print round(100*float(k)/(num_iter), 3)

#time = np.linspace(0, (num_iter-1)*win_shift, num_iter)/samp_rate
#plt.subplot(313)
plt.plot(time, rate, label = 'From code')

# plotting ground truth
offset = 20

#fhr-03_mono
# time_gt = np.asarray([210.15, 210.57, 210.98, 211.40, 211.80, 212.20, 212.61, 213.01, 213.43, 213.83, 214.26, 214.66, 215.04, 215.46, 215.88, 216.27, 216.67, 217.08, 217.49, 217.87, 218.27, 218.7, 219.1, 219.5, 219.91])
# period_gt = np.asarray([416, 415, 421, 394, 403, 412, 402, 414, 403, 429, 402, 384, 412, 426, 386, 397, 411, 416, 375, 405, 429, 404, 399, 404, 406])
# time_gt1 = [227.68, 228.11, 228.52, 228.94, 229.33, 229.67, 230.14, 230.34, 230.55, 230.73, 230.96, 231.34, 231.79, 232.19]
# period_gt1 = [430, 410, 414, 396, 334, 478, 198, 205, 182, 231, 377, 458, 397, 416]
# time_gt2 = [263.02, 263.46, 263.9, 264.3, 264.47, 264.7, 265.12, 265.56, 265.99, 266.4, 266.84, 267.24, 267.67, 268.1, 268.47, 268.9, 269.3, 269.72, 270.13, 270.52, 270.94, 271.34, 271.74, 272.12, 272.55, 272.92, 273.35, 273.7, 274.13, 274.53, 274.92, 275.33, 275.71, 276.1, 276.49, 276.88, 277.28, 277.67, 278.06, 278.45, 278.82, 279.23, 279.63, 279.99, 280.39, 280.75, 281.15]
# period_gt2 = [446, 438, 399, 169, 237, 416, 442, 428, 409, 439, 399, 434, 425, 377, 326, 397, 428, 402, 392, 421, 399, 398, 381, 429, 376, 429, 352, 426, 401, 394, 406, 382, 386, 396, 389, 394, 390, 390, 393, 369, 405, 400, 367, 395, 362, 403, 398]
# time_gt3 = [313.72, 314.1, 314.5, 314.92, 315.34, 315.73, 316.16, 316.55, 316.93, 317.34, 317.75, 318.17, 318.56, 318.95, 319.36, 319.75, 320.14, 320.53, 320.92, 321.29, 321.7]
# period_gt3 = [377, 402, 421, 413, 398, 426, 396, 373, 414, 412, 421, 384, 392, 404, 392, 391, 387, 399, 364, 416, 406]
# time_gt4 = [330.48, 330.86, 331.22, 331.61, 331.94, 332.31, 332.72, 333.08, 333.46, 333.84, 334.21, 334.57, 334.94]
# period_gt4 = [384, 358, 293, 328, 375, 407, 356, 382, 379, 375, 356, 368, 376]
#
# rate_gt = 60.0*1000/period_gt - offset
# rate_gt1 = 60.0*1000/np.asarray(period_gt1) - offset
# rate_gt2 = 60.0*1000/np.asarray(period_gt2) - offset
# rate_gt3 = 60.0*1000/np.asarray(period_gt3) - offset
# rate_gt4 = 60.0*1000/np.asarray(period_gt4) - offset
#
# plt.plot(time_gt, rate_gt, 'r', label = 'Ground truth')
# plt.plot(time_gt1, rate_gt1, 'r')
# plt.plot(time_gt2, rate_gt2, 'r')
# plt.plot(time_gt3, rate_gt3, 'r')
# plt.plot(time_gt4, rate_gt4, 'r')

#fhr-02_mono
# time_gt = [35.04, 35.47, 35.86, 36.3, 36.71, 37.13, 37.54, 37.95, 38.35, 38.77, 39.16, 39.57, 39.98, 40.39, 40.81, 41.21, 41.6, 42.02, 42.42, 42.83, 43.23, 43.64, 44.05, 44.47, 44.88]
# period_gt = [432, 383, 443, 416, 416, 409, 410, 405, 411, 394, 409, 410, 414, 414, 400, 391, 419, 398, 417, 396, 413, 408, 422, 411, 429]
# time_gt1 = [78.69, 79.1, 79.5, 79.92, 80.32, 80.72, 81.13, 81.55, 81.94, 82.35, 82.74, 83.15, 83.54, 83.96, 84.34, 84.75, 85.14, 85.51, 85.92, 86.32, 86.7, 87.09, 87.45, 87.86, 88.24, 88.62, 89, 89.38, 89.79, 90.15, 90.51, 90.89, 91.26, 91.64, 92, 92.39, 92.75, 93.13, 93.49, 93.89, 94.24, 94.62, 94.99, 95.34, 95.71, 96.08, 96.46, 96.82, 97.18, 97.56, 97.93, 98.31, 98.67, 99.04, 99.43, 99.81, 100.23, 100.59, 100.98, 101.41]
# period_gt1 = [410, 397, 418, 402, 397, 415, 415, 395, 404, 388, 410, 391, 421, 381, 408, 388, 377, 406, 399, 387, 381, 369, 406, 384, 371, 382, 381, 407, 365, 361, 379, 372, 372, 370, 384, 359, 383, 354, 400, 359, 373, 373, 354, 364, 376, 375, 363, 361, 375, 373, 375, 368, 372, 380, 383, 418, 359, 393, 427, 412]
#
# rate_gt = 60.0*1000/np.asarray(period_gt) - offset
# rate_gt1 = 60.0*1000/np.asarray(period_gt1) - offset
#
# plt.plot(time_gt, rate_gt, 'r', label = 'Ground truth')
# plt.plot(time_gt1, rate_gt1, 'r')


#NST_Normal
time_gt = [44.37, 44.79, 45.23, 45.66, 46.1, 46.51, 46.92, 47.33, 47.77, 48.17, 48.6, 49, 49.42, 49.82, 50.23, 50.66, 51.06, 51.47, 51.89, 52.31, 52.73, 53.15, 53.58, 54.02, 54.45, 54.91, 55.35, 55.76, 56.23, 56.64, 57.06, 57.52, 57.93, 58.36, 58.82, 59.25, 59.69]
period_gt = [420, 437, 432, 432, 416, 407, 416, 432, 407, 424, 407, 416, 399, 416, 424, 399, 416, 416, 424, 423, 416, 432, 441, 431, 458, 441, 407, 466, 416, 415, 458, 416, 432, 458, 432, 432, 432]
time_gt1 = [115.6, 116.07, 116.52, 116.98, 117.41, 117.86, 118.29, 118.74, 119.16, 119.58, 120, 120.41, 120.84, 121.25, 121.67, 122.04, 122.45, 122.9, 123.3, 123.72, 124.17, 124.58, 124.99, 125.41, 125.87, 126.27, 126.69, 127.15, 127.59, 128.01, 128.43, 128.88, 129.31, 129.73, 130.17, 130.61,
            131.06, 131.52, 131.97, 132.43, 132.89, 133.36, 133.78, 134.24, 134.7, 135.16, 135.59, 136.04, 136.5, 136.96, 137.4, 137.84, 138.29, 138.73, 139.18, 139.63, 140.08, 140.53, 140.99, 141.41, 141.86, 142.3, 142.76, 143.2, 143.65, 144.08, 144.51, 144.94, 145.37, 145.83, 146.24, 146.67]
period_gt1 = [467, 456, 456, 434, 445, 434, 445, 423, 423, 423, 400, 434, 412, 419, 372, 411, 445, 401, 423, 445, 412, 412, 423, 456, 401, 423, 461, 434, 423, 421, 445, 434, 423, 434, 439, 451,
              461, 457, 457, 463, 470, 418, 460, 452, 468, 432, 442, 459, 459, 441, 445, 447, 446, 442, 457, 448, 449, 462, 424, 443, 442, 458, 438, 457, 432, 424, 434, 431, 452, 411, 439, 432]
time_gt2 = [147.11, 147.53, 147.95, 148.38, 148.8, 149.23, 149.69, 150.08, 150.53, 150.94, 151.38, 151.8, 152.23, 152.67, 153.1, 153.55, 153.97, 154.43, 154.87, 155.35, 155.78, 156.22, 156.65, 157.08, 157.52, 157.95, 158.39, 158.81, 159.24, 159.69, 160.13, 160.57, 161.01, 161.45, 161.9, 162.32, 162.76, 163.18, 163.62]
period_gt2 = [424, 424, 429, 414, 434, 455, 394, 444, 414, 434, 424, 434, 434, 434, 445, 424, 455, 445, 475, 434, 434, 434, 434, 434, 434, 433, 426, 431, 443, 441, 445, 439, 440, 449, 420, 440, 423, 438, 440]
time_gt3 = [164.06, 164.48, 164.92, 165.33, 165.76, 166.19, 166.6, 167.05, 167.47, 167.9, 168.32, 168.76, 169.18, 169.61, 170.06, 170.49, 170.94, 171.37, 171.81, 172.23, 172.68, 173.1, 173.55, 174, 174.45, 174.9, 175.37, 175.82, 176.27, 176.71, 177.17, 177.62, 178.08, 178.55, 178.98, 179.44, 179.91, 180.35, 180.82, 181.26, 181.7, 182.15, 182.59, 183.03, 183.47, 183.92, 184.34, 184.79, 185.24, 185.66, 186.11, 186.54, 186.97, 187.41, 187.84, 188.27,
            188.7, 189.16, 189.6, 190.04, 190.49, 190.94, 191.39, 191.81, 192.26, 192.71, 193.13, 193.59, 194.02, 194.46, 194.89, 195.33, 195.76, 196.2, 196.62, 197.08, 197.49, 197.93, 198.35, 198.78, 199.21, 199.67, 200.11, 200.54, 200.98, 201.42, 201.85, 202.29, 202.73, 203.17, 203.6, 204.03, 204.39, 204.91, 205.33, 205.78]
period_gt3 = [416, 443, 406, 433, 433, 404, 449, 424, 429, 422, 438, 422, 431, 450, 428, 449, 428, 439, 424, 449, 424, 445, 452, 449, 452, 472, 446, 449, 449, 458, 446, 460, 468, 429, 458, 470, 441, 472, 438, 436, 449, 449, 435, 438, 449, 424, 449, 449, 428, 443, 433, 432, 433, 433, 428, 438,
              452, 439, 444, 452, 445, 451, 426, 448, 445, 428, 455, 432, 434, 434, 438, 430, 438, 425, 454, 411, 438, 428, 430, 423, 458, 441, 431, 443, 439, 434, 437, 440, 441, 433, 429, 361, 513, 420, 453, 421]
time_gt4 = [0.16, 0.59, 1.02, 1.44, 1.87, 2.3, 2.73, 3.16, 3.59, 4.03, 4.48, 4.93, 5.37, 5.81, 6.26, 6.7, 7.13, 7.58, 8.01, 8.45, 8.89, 9.32, 9.75, 10.18, 10.6, 11.03, 11.45, 11.89, 12.32, 12.79, 13.18, 13.61, 14.01, 14.5, 14.94, 15.36, 15.76, 16.16, 16.61, 17.13, 17.58, 18.02, 18.45, 18.91, 19.31, 19.77, 20.15, 20.58, 21.05, 21.44, 21.84, 22.25, 22.67, 23.14, 23.56, 23.96, 24.38, 24.81, 25.23, 25.67, 26.1, 26.53,
            26.95, 27.36, 27.77, 28.21, 28.63, 29.07, 29.51, 29.92, 30.37, 30.85, 31.27, 31.72, 32.15, 32.57, 33.01, 33.46, 33.91, 34.36, 34.77, 35.23, 35.67, 36.1, 36.53, 36.98, 37.4, 37.84, 38.27, 38.71, 39.13, 39.57, 40, 40.43, 40.84, 41.3, 41.74, 42.17, 42.62, 43.06, 43.49, 43.94]
period_gt4 = [436, 426, 427, 427, 429, 431, 430, 433, 435, 446, 450, 444, 442, 444, 445, 432, 445, 435, 434, 443, 425, 433, 431, 425, 425, 422, 435, 433, 476, 389, 429, 401, 486, 439, 420, 403, 400, 448, 526, 449, 432, 432, 458, 403, 458, 386, 431, 470, 390, 398, 408, 420, 468, 432, 390, 415, 432, 424, 440, 424, 428, 422,
              412, 408, 441, 424, 432, 441, 414, 448, 477, 421, 455, 434, 418, 436, 456, 443, 448, 418, 460, 438, 433, 431, 445, 419, 441, 431, 444, 421, 432, 430, 433, 415, 454, 440, 433, 448, 439, 438, 445, 436]

rate_gt = 60.0*1000/np.asarray(period_gt) - offset
rate_gt1 = 60.0*1000/np.asarray(period_gt1) - offset
rate_gt2 = 60.0*1000/np.asarray(period_gt2) - offset
rate_gt3 = 60.0*1000/np.asarray(period_gt3) - offset
rate_gt4 = 60.0*1000/np.asarray(period_gt4) - offset

plt.plot(time_gt, rate_gt, 'r', label = 'Ground truth')
plt.plot(time_gt1, rate_gt1, 'r')
plt.plot(time_gt2, rate_gt2, 'r')
plt.plot(time_gt3, rate_gt3, 'r')
plt.plot(time_gt4, rate_gt4, 'r')

plt.legend()
plt.ylim([50,250])
plt.grid(True)
plt.xticks(np.arange(tmin, tmin+210+30, 30))
plt.xlabel('Time (in s)')
plt.ylabel('Heart rate (in bpm)')
plt.show()
