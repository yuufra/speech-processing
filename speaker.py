from scipy.io.wavfile import read
import numpy as np
import matplotlib.pyplot as plt
import itertools

speakers = ['spk1','spk2','spk3','spk4']

for s in speakers:
    print('speaker:'+s)
    # average[s] = 0

    for t in itertools.permutations(['a','i','u','e','o']):
        for i in range(1,5):
            # ファイル取得
            print('content: '+''.join(t)+', num: '+str(i))
            wavfile = "./new_vowels/train/"+s+"/"+s+"_"+str(''.join(t))+"_0"+str(i)+".wav"
            # wavfile = "speech_sample/A_a.wav"
            _, data = read(wavfile)

            # 無音区間の除去
            start = 0
            for i in range(len(data)):
                if sum(data[i:i+10]) > 1500:
                    start = i
                    break

            end = len(data)-1
            for i in range(len(data),0,-1):
                if sum(data[i:i+10]) > 1000:
                    end = i+10
                    break

            data = data[start:end]

            # 平均ケプストラムベクトルの算出
            spec = np.fft.fft(data)
            spec = np.log10(abs(spec))

            ceps = np.fft.ifft(spec)
            ceps = ceps[1:]

            # 描画の確認
            # x = [i for i in range(len(ceps))]
            # plt.plot(x,np.abs(ceps))
            # plt.show()