from scipy.optimize import curve_fit
import numpy as np
from scipy.sparse import data
import math
import sys

#validation fail threshold
ERROR_THRESH = 0.025
MEAN_ERROR_THRESH = 0.05

num_voices =6
def func(X, a, b, c,d,e):
    return a  + b*X + c*np.exp2(d*X+e)

arr_0 =[-0.70145603404, 0.0077850969565, 0.0008244905800468867, 11.259986155398293, -0.4975095413274174, 300.0]
def func_simple(X, a, b):
    #X = X/2 -0.5
    global arr_0
    #b = arr_0[1]
    c = arr_0[2]
    d = arr_0[3]
    e = arr_0[4]
    return a  + b*X + c*np.exp2(d*X+e)
sumg = 0
all_popt = []

all_x = []
all_y = []
path_samples = "/udata/delia/tuning/"
path_model= "/udata/delia/calibration/"
step_thresh = 0.0000001
simple_model = False
validate = False
# arguments
n = len(sys.argv)
print("start")
print("Total arguments passed:", n)
if(n > 1):
    val = sys.argv[1]
    print(val)

    if("simple" in sys.argv[1]):
        print("\n run simple analysis")
        simple_model = True
    if("validate" in sys.argv[1]):
        print("\n run validation analysis")
        validate = True
        
exit_code = 0 
# validation
if(validate):
    final_errors = []
    mean_errors = []
    for voicen in range(6):
        filestring  = path_samples + "voice_" + str(voicen) + "_filter.dat"
        audiostring  = path_samples + "voice_" + str(voicen) + "_signal.dat"
        print( filestring)
        stims = ([])
        temps = []
        resf = []
        base_temp = -100
        array = np.fromfile(filestring, dtype="<f")
        audio = np.fromfile(audiostring, dtype="<f")
        temp = array[-1]
        res_array = array[1::3]
        freq_array = array[0::3]
        stim_array = array[2::3]
        
        jumps = np.diff(stim_array)
        jumps = np.absolute(jumps) > 0.0001
        starts = np.array((np.where(jumps==True)))
        stim_slice = stim_array[starts]
        stim_slice = stim_slice.flatten()
        stim_slice = np.append(stim_slice,stim_array[-1])
        starts = starts.flatten()
        starts = np.append(starts,len(freq_array))
        starts = np.insert(starts,0,0)
        audio_starts = starts * len(audio) / len(freq_array)
        start = 0
        sections = len(starts)
        offset = 500
        sample_rate = 1/96000
        f_max_ar = []
        for sec in range(sections-1):
            a_d = audio[math.floor(audio_starts[sec]) + offset:math.floor(audio_starts[sec+1])]
            res = np.fft.fft(a_d)
            freqs = np.fft.fftfreq(res.shape[0],sample_rate)
            q = np.absolute(res)
            q = q[freqs[:] > 0]
            freqs = freqs[freqs[:] > 0]
            i = np.argmax(q)
            #set values which have low amplitude to zero to avoid picking freq of noise
            if(q[i] < 100):
                freqs[i] = 0
            f_max_ar = np.append(f_max_ar,freqs[i])
                
        stim_values = freq_array[starts[1:-1]]
        stim_values  = np.append(stim_values,freq_array[-1])
        stim_values = stim_values[(f_max_ar[:] > 15) & ( f_max_ar[:] < 19000)]
        stim_slice = stim_slice[(f_max_ar[:] > 15) & ( f_max_ar[:] < 19000)]
        f_max_ar = f_max_ar[(f_max_ar[:] > 15) & (f_max_ar[:] < 19000)]
        print("F array:")
        print(f_max_ar)
        
        #middle freq, input of zero should give an FC of this value
        freq_middle = 800
        freq_max = 25600
        freq_log = (np.log2(f_max_ar) -np.log2(freq_middle))/ (np.log2(freq_max) - np.log2(freq_middle))
        freq_log = (freq_log / 2) + 0.5
        
        #get index of val closest to freq middle (800hz)
        index = np.argmin(np.abs(np.array(f_max_ar)-freq_middle))
        
        #find the mean error of the voice
        log_middle = (np.log2(f_max_ar[index]) -np.log2(freq_middle))/ (np.log2(freq_max) - np.log2(freq_middle))
        log_middle = (log_middle/2) + 0.5
        mean_error = log_middle - stim_slice[index]
        
        
        #find the curve error, ignoring the center (mean) error 
        error = freq_log - stim_slice - mean_error
        error = error[(stim_slice[:] > 0.25) & (stim_slice[:] < 1) & (f_max_ar[:] > 150)]
        error = np.sqrt(np.mean(error**2))
        final_errors = np.append(final_errors, error) 
        mean_errors = np.append(mean_errors, mean_error)
        
    print("MSE error")
    print(final_errors)
    print("Middle error")
    mean_errors = abs(mean_errors)
    print(mean_errors)
    max_error = np.max(final_errors)
    max_error_mean = np.max(mean_errors)
    exit_code = 0
    if(max_error > ERROR_THRESH):
        exit_code = 1
        print("VALIDATE FAILED")
    elif (max_error_mean > MEAN_ERROR_THRESH):
        exit_code = 1
        print("VALIDATE FAILED")
    else:
        print("PASS")
    
    with open(path_model + "filter_cal_status.txt", 'w') as f:
        f.write(str(exit_code) + '\n')
        for i in range(len(final_errors)):
            f.write("{:5.4f}".format(final_errors[i]) + ' ')
else:
    for voicen in range(6):
        filestring   = path_samples + "voice_" + str(voicen) + "_filter.dat"
        audiostring  = path_samples + "voice_" + str(voicen) + "_signal.dat"
        print( filestring)
        stims = ([])
        temps = []
        resf = []
        base_temp = -100
        array = np.fromfile(filestring, dtype="<f")
        audio = np.fromfile(audiostring, dtype="<f")
        temp = array[-1]
        print("temp: "+str(temp))
        array = array[0:-1]
        res_array = array[1::3]
        freq_array = array[0::3]
        stim_array = array[2::3]
        jumps = np.diff(stim_array)
        jumps = np.absolute(jumps) > step_thresh
        starts = np.array((np.where(jumps==True)))
        starts = starts.flatten()
        starts = np.append(starts,len(freq_array))
        starts = np.insert(starts,0,0)
        audio_starts = starts * len(audio) / len(freq_array)
        start = 0
        sections = len(starts)
        print(sections)
        offset = 500
        sample_rate = 1/96000
        f_max_ar = []
        for sec in range(sections-1):
            a_d = audio[math.floor(audio_starts[sec]) + offset:math.floor(audio_starts[sec+1])]
            res = np.fft.fft(a_d)
            freqs = np.fft.fftfreq(res.shape[0],sample_rate)
            q = np.absolute(res)
            q = q[freqs[:] > 0]
            freqs = freqs[freqs[:] > 0]
            i = np.argmax(q)
            f_max_ar = np.append(f_max_ar,freqs[i])
        stim_values = freq_array[starts[1:-1]]
        stim_values  = np.append(stim_values,freq_array[-1])
        stim_values = stim_values[(f_max_ar[:] > 15) & ( f_max_ar[:] < 19000)]
        f_max_ar = f_max_ar[(f_max_ar[:] > 15) & (f_max_ar[:] < 19000)]
        print("F array:")
        print(f_max_ar)
        sigma_in = np.ones(np.shape(stim_values))
        
        popt0  = 1,1,1,1,1
        if(simple_model):
            popt0 = 1,1
    
        #linear model max freq, in practice, the filter will not reach this value 
        freq_max = 25600
        
        #middle freq, input of zero should give an FC of this value
        freq_middle = 800
        freq_log = (np.log2(f_max_ar) -np.log2(freq_middle))/ (np.log2(freq_max) - np.log2(freq_middle))
        freq_log = (freq_log / 2) + 0.5
        all_x  = np.append(all_x, stim_values)
        if(simple_model):
            popt1, pcov = curve_fit(func_simple, freq_log, stim_values,p0 = popt0,sigma=sigma_in ,maxfev=1000000)
        else:
            popt1, pcov = curve_fit(func, freq_log, stim_values,p0 = popt0,sigma=sigma_in ,maxfev=1000000)
        print("popt voice  " + str(voicen))
        print(popt1)
        all_y = np.append(all_y, freq_log - popt1[0])
        file_items = 9
        model_arr  = np.zeros(file_items)
        if(simple_model):
            model_arr[0] = popt1[0] 
            model_arr[1] = popt1[1] 
            model_arr[2] = arr_0[2]
            model_arr[3] = arr_0[3]
            model_arr[4] = arr_0[4]
        else:
            model_arr[0] = popt1[0]
            model_arr[1] = popt1[1]
            model_arr[2] = popt1[2]
            model_arr[3] = popt1[3]
            model_arr[4] = popt1[4]
        with open(path_model + "voice_" +str(voicen) +"_filter.model", 'w') as f:
            for i in range(file_items):
                f.write("{:10.10f}".format(model_arr[i]) + ' ')
            f.write("{:10.10f}".format(base_temp)  + ' ')
        if(simple_model):
            y_est = func_simple((freq_log), *popt1)
        else:
            y_est = func((freq_log), *popt1)
        error = y_est - stim_values

    all_popt, pcov = curve_fit(func, all_y, all_x ,maxfev=1000000)
    print("all popt")
    print(all_popt)
print("done")
sys.exit(exit_code)