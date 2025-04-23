from astropy.io import ascii
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math

def SED_fit(freq, flux, idx):
    freq_log = []
    flux_log = []
    for k in range(len(flux)):
        freq_t = freq[k][idx]
        flux_t = flux[k][idx].ravel()
        zero_idx=np.where(flux_t<0)[0]
        freq_t = np.delete(freq_t,zero_idx)
        flux_t = np.delete(flux_t,zero_idx)
        freq_log.extend(np.log(freq_t).tolist())
        flux_log.extend(np.log(flux_t).tolist())
    
    z, cov = np.polyfit(freq_log, flux_log, 1, cov=True)
    p = np.poly1d(z)
    xp = np.linspace(np.log(freq[0][idx][0]), np.log(freq[4][idx][-1]), 10)
    ax.plot(np.exp(xp), np.exp(p(xp)),'-', color='grey')
    return round(z[0],2)

tracks = ['track1', 'track2']
colors = cm.jet([0.2, 0.3, 0.5, 0.7, 0.8, 0.95])
freq1 = [199, 219, 228, 248]
freq2 = [262, 278, 294, 310]
#freq3 = [337, 357, 399.5, 415.5]
freq = []
flux = []
rms = []
flux_sel = []
rms_sel = []
selcal_target = []
others=['AS_206', 'DoAr_16', 'DoAr_25', 'DoAr_33', 'DoAr_44', 'GSS_26', 'GSS_39', 'HBC_266', 'IRS_37', 'IRS_39', 'IRS_41', 'IRS_51', 'VSSG_1', 'WSB_31', 'WSB_60', 'YLW_8', 'YLW_47']
#others=['AS_206', 'DoAr_16', 'DoAr_33', 'HBC_266', 'IRS_37', 'IRS_39', 'IRS_41', 'IRS_51', 'VSSG_1', 'YLW_47']
others=['DoAr_25', 'DoAr_44', 'GSS_26', 'GSS_39', 'WSB_31', 'WSB_60', 'YLW_8']
# all unresolved sources
Ophidict1 = {
    "DoAr 16": "",
    "YLW 47": "ODISEA_C4_100",
    "IRS 39": "ODISEA_C4_084",
    "DoAr 33": "ODISEA_C4_103",
    "IRS 37": "ODISEA_C4_082",
    "IRS 41": "ODISEA_C4_085",
    "HBC 266": "ODISEA_C4_117",
    "AS 206": "ODISEA_C4_027",
    "WSB 60": "ODISEA_C4_114",
    "DoAr 24E": "ODISEA_C4_037",
    "IRS 51": "ODISEA_C4_105A",
    "DoAr 44": "ODISEA_C4_127",
    "VSSG 1": "ODISEA_C4_034",
    "YLW 8": "ODISEA_C4_027",
    "GSS 26": "ODISEA_C4_030",
    "DoAr 25": "ODISEA_C4_039",
    "GSS 39": "ODISEA_C4_051",
    "WSB 31": "ODISEA_C4_041"
}

Ophidict2 = {
    "DoAr 16": "DoAr_16",
    "YLW 47": "YLW_47" ,
    "IRS 39": "IRS 39, WL 4",
    "DoAr 33": "DoAr_33",
    "IRS 37": "IRS 37, GY 244",
    "IRS 41": "IRS 41, WL 3",
    "HBC 266": "SR_13",
    "AS 206": "SR_4",
    "WSB 60": "WSB_60",
    "DoAr 24E": "DoAr_24E",
    "IRS 51": "IRS 51, GY 315",
    "DoAr 44": "DoAr_44",
    "VSSG 1": "VSSG_1",
    "YLW 8": "SR_21",
    "GSS 26": "GSS_26",
    "DoAr 25": "DoAr_25",
    "GSS 39": "EL_27",
    "WSB 31": "EL_24"
}


Ophidict = {
    "DoAr 16": "DoAr 16, HBC 257",
    "YLW 47": "IRS 49, GY 308",
    "IRS 39": "IRS 39, WL 4",
    "DoAr 33": "DoAr 33, WSB 53",
    "IRS 37": "IRS 37, GY 244",
    "IRS 41": "IRS 41, WL 3",
    "HBC 266": "SR 13, HBC 266",
    "AS 206": "SR 4, AS 206",
    "WSB 60": "WSB 60, YLW 58",
    "DoAr 24E": "DoAr 24E, GSS 31",
    "IRS 51": "IRS 51, GY 315",
    "DoAr 44": "DoAr 44, HBC 268",
    "VSSG 1": "EL 20, VSSG 1",
    "YLW 8": "SR 21, EL 30",
    "GSS 26": "GSS 26",
    "DoAr 25": "DoAr 25, WSB 29",
    "GSS 39": "EL 27, GSS 39",
    "WSB 31": "EL 24, WSB 31"
}


target = []
spidx = []
freq_q = []
Flux_q = []
Flux_e_q = []

freq_l_q = []
Flux_l_e_q = []
Flux_l_q = []
table1 = ascii.read('datafile1.txt')
table2 = np.loadtxt('datafile2.txt',dtype=str)
table3 = np.loadtxt('datafile3.txt',dtype=str)
table4 = np.loadtxt('datafile4.txt',dtype=str)

freq2_q = []
Flux2_q = []
Flux2_e_q = []

freq3_q = []
Flux3_q = []
Flux3_e_q = []

freq4_q = []
Flux4_q = []
Flux4_e_q = []


filename='Fnu_'+tracks[0]+'.txt'
file = open(filename, 'r')
lines = file.readlines()
for i in range(len(lines)):
    if lines[i].split()[0] in others:
        target.append(lines[i].split()[0])
        idx1=np.where(table1['Name'] == Ophidict[lines[i].split()[0].replace('_',' ')])[0]
        idx2=np.where(table2[:,0] == Ophidict2[lines[i].split()[0].replace('_',' ')])[0]
        idx3=np.where(table3[:,0] == Ophidict1[lines[i].split()[0].replace('_',' ')])[0]
        idx4=np.where(table4[:,1] == Ophidict2[lines[i].split()[0].replace('_',' ')])[0]

        print(idx4)

        tempf = []
        tempF = []
        tempF_e = []
        tempf_l = []
        tempF_l = []
        tempF_l_e = []

        if len(idx1) != 0:
            spidx.append(table1['Sp+Index'][idx1][0])
            if table1['l_F850'][idx1][0] == '<':
                tempf_l.append(353)
                tempF_l.append(table1['F850'][idx1][0])
                tempF_l_e.append(table1['e_F850'][idx1][0])
            else:
                tempf.append(353)
                tempF.append(table1['F850'][idx1][0])
                tempF_e.append(table1['e_F850'][idx1][0])

            if table1['l_F1.3'][idx1][0] == '<':
                tempf_l.append(231)
                tempF_l.append(table1['F1.3'][idx1][0])
                tempF_l_e.append(table1['e_F1.3'][idx1][0])
            else:
                tempf.append(231)
                tempF.append(table1['F1.3'][idx1][0])
                tempF_e.append(table1['e_F1.3'][idx1][0])

        else:
            spidx.append(np.nan)
                    
        freq_q.append(tempf)
        Flux_q.append(tempF)
        Flux_e_q.append(tempF_e)
        freq_l_q.append(tempf_l)
        Flux_l_q.append(tempF_l)
        Flux_l_e_q.append(tempF_l_e)

        tempf = []
        tempF = []
        tempF_e = []

        if len(idx2) != 0:
            tempf.append(float(345))
            tempF.append(float(table2[idx2,1][0]))
            tempF_e.append(float(table2[idx2,3][0]))

        freq2_q.append(tempf)
        Flux2_q.append(tempF)
        Flux2_e_q.append(tempF_e)
       
        tempf = []
        tempF = []
        tempF_e = []

        if len(idx3) != 0:
            tempf.append(float(233))
            tempF.append(float(table3[idx3,3]))
            tempF_e.append(float(table3[idx3,4]))

        freq3_q.append(tempf)
        Flux3_q.append(tempF)
        Flux3_e_q.append(tempF_e)

        tempf = []
        tempF = []
        tempF_e = []

        if len(idx4) != 0:
            tempf.append(float(348.5))
            tempF.append(float(table4[idx4,17]))
            tempF_e.append(float(table4[idx4,19]))

        freq4_q.append(tempf)
        Flux4_q.append(tempF)
        Flux4_e_q.append(tempF_e)



for track in tracks:
    filename='freq_'+track+'.txt'
    file = open(filename, 'r')
    lines = file.readlines()
    tempf = []
    a = 0
    for i in range(len(target)):
        for j in range(len(lines)):
            if lines[j].split()[0] == target[i]:
                a = 1
                tempf.append(lines[j].split()[1:])
                break
    if a == 1:
        freq.append(tempf)


for track in tracks:
    filename='Fnu_'+track+'.txt'
    file = open(filename, 'r')
    lines = file.readlines()
    temp = []
#    tempf = []
    a = 0
    for i in range(len(target)):
        for j in range(len(lines)):
            if lines[j].split()[0] == target[i]:
                a = 1
                temp.append(lines[j].split()[1:])
#                if track =='track1' or track == 'track2' or track == 'track12':
#                    tempf.append(freq1)
#                elif track == 'track3':
#                    tempf.append(freq2)
#                else:
#                    tempf.append(freq3)
                break
    if a == 1:
        flux.append(temp)
        freq.append(tempf)

for track in tracks:
    filename='Fnu_e_'+track+'.txt'
    file = open(filename, 'r')
    lines = file.readlines()
    temp = []
    a = 0
    for i in range(len(target)):
        for j in range(len(lines)):
            if lines[j].split()[0] == target[i]:
                a = 1
                temp.append(lines[j].split()[1:])
                break
    if a == 1:
        rms.append(temp)

#filename='../../spectrum_fit/spectrum_fitting_200_420_com/fitting_result/alpha.txt'
#file = open(filename, 'r')
#lines = file.readlines()
#for i in range(len(target)):
#    for j in range(len(lines)):
#        if lines[j].split()[0] == target[i]:
#            alpha.append(float(lines[j].split()[1]))
#            alpha_p.append(float(lines[j].split()[2]))
#            alpha_m.append(-float(lines[j].split()[3]))
#            break

freq=np.array(freq).astype(float)*10**(-9)
flux=np.array(flux).astype(float) 
flux=np.where(flux==0.0,-100,flux)
rms=np.array(rms).astype(float)
#flux_sel=np.array(flux_sel).astype(float)
#rms_sel=np.array(rms_sel).astype(float)*1000

def get_idx(target, flux_0, rank):
    # index [:,0] for 345 GHz lsb
    flux_0_s = sorted(flux_0[:,0], reverse = True)
    idx_l= list(np.where(flux_0[:,0]==flux_0_s[rank])[0])
    return idx_l

num = 9
num_sel = 3

num_figure= int((len(target))/num)

times = 1
for i in range(num_figure):
    fig = plt.figure(figsize=(20, 18))
    for j in range(num):
        ax = fig.add_subplot(331+j)
        rank = num*i+j
        idx_l = get_idx(target, flux[tracks.index('track1')], rank)
        
        if len(idx_l) > times:
            idx = idx_l[times-1]
            times += 1
        elif len(idx_l) == times and times > 1:
            idx = idx_l[times-1]
            times = 1
        else:
            idx = idx_l[0]
            times = 1
        if_zero = 0
        flux_max = 0
#        for k in range(6):
#            ax.errorbar(freq[k][idx], flux[k][idx], yerr=rms[k][idx], fmt='o', color=colors[k], ms=6,alpha=0.8)
#            if_zero+=flux[k][idx].tolist().count(-100)
#        flux_max = max(flux_max, max(flux[:][idx]))

#        freq_com=np.concatenate((freq[0][idx][0:4],freq[1][idx][0:4],freq[2][idx][0:4],freq[3][idx][0:4],freq[4][idx][0:4],freq[5][idx][0:4]))
#        flux_com=np.concatenate((flux[0][idx][0:4],flux[1][idx][0:4],flux[2][idx][0:4],flux[3][idx][0:4],flux[4][idx][0:4],flux[5][idx][0:4]))
#        rms_com=np.concatenate((rms[0][idx][0:4],rms[1][idx][0:4],rms[2][idx][0:4],rms[3][idx][0:4],rms[4][idx][0:4],rms[5][idx][0:4]))

        freq_com=np.concatenate((freq[0][idx][0:4],freq[1][idx][0:4]))
        flux_com=np.concatenate((flux[0][idx][0:4],flux[1][idx][0:4]))
        rms_com=np.concatenate((rms[0][idx][0:4],rms[1][idx][0:4]))
#        freq_com=np.concatenate((freq[0][idx][0:4],freq[1][idx][0:4],freq[2][idx][0:4]))
#        flux_com=np.concatenate((flux[0][idx][0:4],flux[1][idx][0:4],flux[2][idx][0:4]))
#        rms_com=np.concatenate((rms[0][idx][0:4],rms[1][idx][0:4],rms[2][idx][0:4]))

        with open(target[idx].replace(" ","_")+'.specdata.txt', 'w') as f:
            for k in range(len(freq_com)):
                f.writelines(str(freq_com[k])+'   '+str(flux_com[k])+'   '+str(rms_com[k])+'\n')

        ax.errorbar(freq_com, flux_com, yerr=rms_com, fmt='o', color='black',ms=8,alpha=0.8)
#        legend=['1126 230 GHz-1','1018 230 GHz-2','1121 270 GHz-3','0925 400 GHz-4','0903 400 GHz-5','1122 400 GHz-6','com_vis']
        legend=['This work']
        flux_max = max(flux_max, max(flux_com))

        if len(Flux_q[idx]) > 0:
            try:
                ax.errorbar(freq_q[idx], Flux_q[idx], yerr=Flux_e_q[idx], fmt='D', color='grey',alpha=0.5)
                flux_max = max(flux_max, max(Flux_q[idx]))
            except:
                print(freq_q[idx], Flux_q[idx],Flux_e_q[idx])
            legend.append('A&W(2007)')

        if len(Flux_l_q[idx]) > 0:
            ax.errorbar(freq_q[idx], Flux_q[idx], yerr=Flux_e_q[idx], fmt='D', color='grey',alpha=0.5)
            flux_max = max(flux_max, max(Flux_l_q[idx]))
            legend.append('A&W(2007)')
        
        if len(Flux2_q[idx]) > 0:
            try:
                ax.errorbar(freq2_q[idx], Flux2_q[idx], yerr=Flux2_e_q[idx], fmt='D', color='navy',alpha=0.5)
                flux_max = max(flux_max, max(Flux2_q[idx]))
            except:
                print(freq2_q[idx], Flux2_q[idx],Flux2_e_q[idx])
            legend.append('Andrews(2009,2010)')

        if len(Flux3_q[idx]) > 0:
            try:
                ax.errorbar(freq3_q[idx], Flux3_q[idx], yerr=Flux3_e_q[idx], fmt='D', color='turquoise',alpha=0.5)
                flux_max = max(flux_max, max(Flux3_q[idx]))
            except:
                print(freq3_q[idx], Flux3_q[idx],Flux3_e[idx])
            legend.append('Cieza(2019)')
#            with open(target[idx].replace(" ","_")+'.specdata.txt', 'a') as f:
#                for k in range(len(Flux3_q[idx])):
#                    f.writelines(str(freq3_q[idx][k])+'   '+str(Flux3_q[idx][k])+'   '+str(Flux3_e_q[idx][k])+'\n')

        if len(Flux4_q[idx]) > 0:
            try:
                ax.errorbar(freq4_q[idx], Flux4_q[idx], yerr=Flux4_e_q[idx], fmt='D', color='limegreen',alpha=0.5)
                flux_max = max(flux_max, max(Flux4_q[idx]))
            except:
                print(freq4_q[idx], Flux4_q[idx],Flux4_e[idx])
            legend.append('Cox(2017)')
#            with open(target[idx].replace(" ","_")+'.specdata.txt', 'a') as f:
#                for k in range(len(Flux4_q[idx])):
#                    f.writelines(str(freq4_q[idx][k])+'   '+str(Flux4_q[idx][k])+'   '+str(Flux4_e_q[idx][k])+'\n')


        plt.ylim([0, flux_max*1.25])
        ax.tick_params(axis='both', which='major', labelsize=12)
        plt.legend(legend,ncol=2,loc=2,fontsize=10)

        
        if (j%3==0):
            plt.ylabel('$F_\\nu$  (mJy)', size=14)
        if (j/3>0):
            plt.xlabel('$\\nu$  (GHz)', size=14) 
#        if (if_zero>0):
#            ax.title.set_text(target[idx])
#        else:
#        alpha = SED_fit(freq, flux, idx)
#        ax.title.set_text('%s $\u03B1 = %.2f \pm ^{%.2f} _{%.2f}$ , [%.2f]' % (target[idx].replace('_',' '), alpha[idx], alpha_p[idx], alpha_m[idx], spidx[idx]))

#        ax.title.set_text('%s' % (target[idx].replace('_',' ')))
        ax.title.set_text(target[idx].replace('_',' '))
        ax.title.set_size(20)
    fig.tight_layout()
    plt.savefig('flux_measurement_0_'+str(i)+'.pdf', format='PDF', transparent=True)
    plt.savefig('flux_measurement_0_'+str(i)+'.png', transparent=True)
    plt.close(fig) 



i = num_figure
fig = plt.figure(figsize=(20, 18))
for j in range(len(target)%num):

        ax = fig.add_subplot(331+j)
        rank = num*i+j
        idx_l = get_idx(target, flux[tracks.index('track1')], rank)

        if len(idx_l) > times:
            idx = idx_l[times-1]
            times += 1
        elif len(idx_l) == times and times > 1:
            idx = idx_l[times-1]
            times = 1
        else:
            idx = idx_l[0]
            times = 1
        if_zero = 0
        flux_max = 0
#        for k in range(6):
#            ax.errorbar(freq[k][idx], flux[k][idx], yerr=rms[k][idx], fmt='o', color=colors[k], ms=6,alpha=0.8)
#            if_zero+=flux[k][idx].tolist().count(-100)
#        flux_max = max(flux_max, max(flux[:][idx]))

#        freq_com=np.concatenate((freq[0][idx][0:4],freq[1][idx][0:4],freq[2][idx][0:4],freq[3][idx][0:4],freq[4][idx][0:4],freq[5][idx][0:4]))
#        flux_com=np.concatenate((flux[0][idx][0:4],flux[1][idx][0:4],flux[2][idx][0:4],flux[3][idx][0:4],flux[4][idx][0:4],flux[5][idx][0:4]))
#        rms_com=np.concatenate((rms[0][idx][0:4],rms[1][idx][0:4],rms[2][idx][0:4],rms[3][idx][0:4],rms[4][idx][0:4],rms[5][idx][0:4]))

#        freq_com=np.concatenate((freq[0][idx][0:4],freq[1][idx][0:4],freq[2][idx][0:4]))
#        flux_com=np.concatenate((flux[0][idx][0:4],flux[1][idx][0:4],flux[2][idx][0:4]))
#        rms_com=np.concatenate((rms[0][idx][0:4],rms[1][idx][0:4],rms[2][idx][0:4]))
#        freq_com=freq[0][idx][0:4]
#        flux_com=flux[0][idx][0:4]
#        rms_com=rms[0][idx][0:4]
        freq_com=np.concatenate((freq[0][idx][0:4],freq[1][idx][0:4]))
        flux_com=np.concatenate((flux[0][idx][0:4],flux[1][idx][0:4]))
        rms_com=np.concatenate((rms[0][idx][0:4],rms[1][idx][0:4]))
        with open(target[idx].replace(" ","_")+'.specdata.txt', 'w') as f:
            for k in range(len(freq_com)):
                f.writelines(str(freq_com[k])+'   '+str(flux_com[k])+'   '+str(rms_com[k])+'\n')

        ax.errorbar(freq_com, flux_com, yerr=rms_com, fmt='o', color='black',ms=8,alpha=0.8)
#        legend=['1126 230 GHz-1','1018 230 GHz-2','1121 270 GHz-3','0925 400 GHz-4','0903 400 GHz-5','1122 400 GHz-6','com_vis']
        legend=['This work']
        flux_max = max(flux_max, max(flux_com))

        if len(Flux_q[idx]) > 0:
            try:
                ax.errorbar(freq_q[idx], Flux_q[idx], yerr=Flux_e_q[idx], fmt='D', color='grey',alpha=0.5)
                flux_max = max(flux_max, max(Flux_q[idx]))
            except:
                print(freq_q[idx], Flux_q[idx],Flux_e_q[idx])
            legend.append('A&W(2007)')

        if len(Flux_l_q[idx]) > 0:
            ax.errorbar(freq_q[idx], Flux_q[idx], yerr=Flux_e_q[idx], fmt='D', color='grey',alpha=0.5)
            flux_max = max(flux_max, max(Flux_l_q[idx]))
            legend.append('A&W(2007)')

        if len(Flux2_q[idx]) > 0:
            try:
                ax.errorbar(freq2_q[idx], Flux2_q[idx], yerr=Flux2_e_q[idx], fmt='D', color='navy',alpha=0.5)
                flux_max = max(flux_max, max(Flux2_q[idx]))
            except:
                print(freq2_q[idx], Flux2_q[idx],Flux2_e_q[idx])
            legend.append('Andrews(2009,2010)')

        if len(Flux3_q[idx]) > 0:
            try:
                ax.errorbar(freq3_q[idx], Flux3_q[idx], yerr=Flux3_e_q[idx], fmt='D', color='turquoise',alpha=0.5)
                flux_max = max(flux_max, max(Flux3_q[idx]))
            except:
                print(freq3_q[idx], Flux3_q[idx],Flux3_e[idx])
            legend.append('Cieza(2019)')
#            with open(target[idx].replace(" ","_")+'.specdata.txt', 'a') as f:
#                for k in range(len(Flux3_q[idx])):
#                    f.writelines(str(freq3_q[idx][k])+'   '+str(Flux3_q[idx][k])+'   '+str(Flux3_e_q[idx][k])+'\n')

        if len(Flux4_q[idx]) > 0:
            try:
                ax.errorbar(freq4_q[idx], Flux4_q[idx], yerr=Flux4_e_q[idx], fmt='D', color='limegreen',alpha=0.5)
                flux_max = max(flux_max, max(Flux4_q[idx]))
            except:
                print(freq4_q[idx], Flux4_q[idx],Flux4_e[idx])
            legend.append('Cox(2017)')
#            with open(target[idx].replace(" ","_")+'.specdata.txt', 'a') as f:
#                for k in range(len(Flux4_q[idx])):
#                    f.writelines(str(freq4_q[idx][k])+'   '+str(Flux4_q[idx][k])+'   '+str(Flux4_e_q[idx][k])+'\n')


        plt.ylim([0, flux_max*1.25])
        ax.tick_params(axis='both', which='major', labelsize=12)
        plt.legend(legend,ncol=2,loc=2,fontsize=10)


        if (j%3==0):
            plt.ylabel('$F_\\nu$  (mJy)', size=14)
        if (j/3>0):
            plt.xlabel('$\\nu$  (GHz)', size=14)
#        if (if_zero>0):
#            ax.title.set_text(target[idx])
#        else:
#        alpha = SED_fit(freq, flux, idx)
#        ax.title.set_text('%s $\u03B1 = %.2f \pm ^{%.2f} _{%.2f}$ , [%.2f]' % (target[idx].replace('_',' '), alpha[idx], alpha_p[idx], alpha_m[idx], spidx[idx]))

#        ax.title.set_text('%s' % (target[idx].replace('_',' ')))
        ax.title.set_text(target[idx].replace('_',' '))
        ax.title.set_size(20)
fig.tight_layout()
plt.savefig('flux_measurement_0_'+str(i)+'.pdf', format='PDF', transparent=True)
plt.savefig('flux_measurement_0_'+str(i)+'.png', transparent=True)
plt.close(fig)

