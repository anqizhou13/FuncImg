def max_value(inputlist):
    return max([sublist[-1] for sublist in inputlist])

def min_value(inputlist):
    return min([sublist[-1] for sublist in inputlist])

def reject_outliers(data, m = 2.):
    import numpy as np
    data = np.array(data)
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return data[s<m]

def process(input,paths):
    import glob
    import os
    import numpy as np
    import pandas as pd
    import pickle

    totalTime = 15

    for path in paths:
        #for each path
        files = glob.glob("{}/{}/**/Values.csv".format(input,path), recursive = True)

        raw = []
        dF = []
        dF_f0 = []
        mean_baseline = []
        mean_dF_f0 = []
        max_dF_f0 = []

        for file in files:

            with open(file) as f:
                lines = (line for line in f if not line.startswith('#'))
                data = np.loadtxt(lines, delimiter=',', skiprows=1)
                data = np.transpose(data)   
                data = data[-1]

                if len(data) > totalTime:
                    data = data[0:totalTime-1]

                raw.append(data)
            
                # mean baseline value is the average of the first 5 elements
                baseline = np.mean(data[0:4]) 
                mean_baseline.append(baseline)
                dF.append(data-baseline)
                val = (data-baseline)/baseline
                dF_f0.append(val)
                mean_dF_f0.append(np.mean(val[5:9]))
                max_dF_f0.append(np.array(val[5:9]).max())

        pickle.dump([raw,dF,dF_f0,mean_baseline,mean_dF_f0, max_dF_f0],open('{}/{}/Values.pkl'.format(input,path),'wb'))

def process_dualChannels(input,paths):
    import glob
    import os
    import numpy as np
    import pandas as pd
    import pickle

    totalTime = 15

    for path in paths:
        #for each path

        ## process data files from channel 1
        output = "{}/{}/C1".format(input,path)
        if not os.path.exists(output):
            os.makedirs(output)

        files = glob.glob("{}/{}/**/Values_C1.csv".format(input,path), recursive = True)

        raw = []
        dF = []
        dF_f0 = []
        mean_baseline = []
        mean_dF_f0 = []
        max_dF_f0 = []
        channel = []
        genotype = []
        mean_baseline_all = []
        mean_dF_F0_all = []
        trial = []

        for file in files:

            with open(file) as f:
                lines = (line for line in f if not line.startswith('#'))
                data = np.loadtxt(lines, delimiter=',', skiprows=1)
                data = np.transpose(data)   
                data = data[-1]

                if len(data) > totalTime:
                    data = data[0:totalTime-1]

                baseline = np.mean(data[-5:]) 
                data = np.concatenate((data[-5:],data[:len(data)-5]),axis=None)
                raw.append(data)
            
                # mean baseline value is the average of the first 5 elements
                mean_baseline.append(baseline)
                dF.append(data-baseline)
                val = (data-baseline)/baseline
                dF_f0.append(val)
                mean_dF_f0.append(np.mean(val[5:9]))
                max_dF_f0.append(np.array(val[5:9]).max())
                channel.append('syp-pHTomato')
                genotype.append(path)
                mean_baseline_all.append(baseline)
                mean_dF_F0_all.append(np.mean(val[5:9]))
                trial.append(len(mean_baseline))

        features_C1 = [raw,dF,dF_f0,mean_baseline,mean_dF_f0, max_dF_f0]
        pickle.dump(features_C1,open('{}/Values.pkl'.format(output),'wb'))

        ## process data files from channel 2
        output = "{}/{}/C2".format(input,path)
        if not os.path.exists(output):
            os.makedirs(output)

        files = glob.glob("{}/{}/**/Values_C2.csv".format(input,path), recursive = True)

        raw = []
        dF = []
        dF_f0 = []
        mean_baseline = []
        mean_dF_f0 = []
        max_dF_f0 = []

        for file in files:

            with open(file) as f:
                lines = (line for line in f if not line.startswith('#'))
                data = np.loadtxt(lines, delimiter=',', skiprows=1)
                data = np.transpose(data)   
                data = data[-1]

                if len(data) > totalTime:
                    data = data[0:totalTime-1]

                baseline = np.mean(data[-5:]) 
                data = np.concatenate((data[-5:],data[:len(data)-5]),axis=None)
                raw.append(data)
            
                # mean baseline value is the average of the first 5 elements
                mean_baseline.append(baseline)
                dF.append(data-baseline)
                val = (data-baseline)/baseline
                dF_f0.append(val)
                mean_dF_f0.append(np.mean(val[5:9]))
                max_dF_f0.append(np.array(val[5:9]).max())
                channel.append('syp-GCaMP3')
                genotype.append(path)
                mean_baseline_all.append(baseline)
                mean_dF_F0_all.append(np.mean(val[5:9]))
                trial.append(len(mean_baseline))


        features_C2 = [raw,dF,dF_f0,mean_baseline,mean_dF_f0, max_dF_f0]
        pickle.dump(features_C2,open('{}/Values.pkl'.format(output),'wb'))

        ## compute ratio between C1 and C2
        df = pd.DataFrame(data = np.transpose(np.array([trial,mean_baseline_all,mean_dF_F0_all,channel,genotype])), columns = ['Trial','Baseline F0','Mean dF/F0','Channel','Genotype'])
        df.to_csv('{}/{}/summary.csv'.format(input,path))



def analyze(input,output,genotypes,colors):
    import pickle
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    from itertools import combinations
    import scipy.stats as stats

    if not os.path.exists(output):
        os.makedirs(output)

    data = []
    baseline_all = []
    dF_F0_all = []
    metrics = ['Raw ROI','dF','dF_F0']

    totalTime = 15

    for genotype in genotypes:
        raw,dF,dF_f0,mean_baseline,mean_dF_f0 = pickle.load(open('{}/{}/Values.pkl'.format(input,genotype), 'rb'))
        
        vars = [raw,dF,dF_f0]

        # for each metric
        for j,var in enumerate(vars):

            mean = []
            sem = []
            # convert the variable to array
            var = np.array(var)

            # then for each time point, compute the mean and sem over all larvae
            for t in range(totalTime-1):
                mean.append(np.mean([var[l][t] for l in range(len(var))]))
                sem.append(stats.sem([var[l][t] for l in range(len(var))]))
            
            # append to the dataframe
            data.append([mean,sem])

        
        baseline_all.append(mean_baseline)
        dF_F0_all.append(mean_dF_f0)
        
    # convert the data to array    
    data = np.array(data)
    # reshape the data
    dims = data.shape
    data = data.reshape(len(genotypes),len(metrics),dims[-2],dims[-1])

    # data is then a 3D list
    # 1st dim: genotype or conditions
    # 2nd dim: always of size 3, raw, dF, and dF_f0
    # 3rd dim: always of size 2, mean and sem

    # visualize
    # Say, "the default sans-serif font is COMIC SANS"
    plt.rcParams['font.sans-serif'] = "Arial"
    # Then, "ALWAYS use sans-serif fonts"
    plt.rcParams['font.family'] = "sans-serif"

    fig,axs = plt.subplots(1,3, figsize=(11,8), dpi=300)
    x = np.linspace(0,totalTime-1,totalTime-1)
    for i,metric in enumerate(metrics):
        for j,__ in enumerate(genotypes):
            axs[i].plot(x,data[j][i][0], colors[j])
            axs[i].fill_between(x,np.array(data[j][i][0])-np.array(data[j][i][1]),np.array(data[j][i][0])+np.array(data[j][i][1]), facecolor = colors[j], alpha = 0.3)
            axs[i].set_title(metric)
            axs[i].set_xlabel('Seconds (s)')

            x_left, x_right = axs[i].get_xlim()
            y_low, y_high = axs[i].get_ylim()       
            axs[i].set_aspect(abs((x_right-x_left)/(y_low-y_high))*1)
            axs[i].spines.right.set_visible(False)
            axs[i].spines.top.set_visible(False)
    plt.savefig('{}/lineplot.svg'.format(output), format='svg')
    plt.savefig('{}/lineplot.pdf'.format(output), format='pdf')


    result = []
    combos = list(combinations(range(len(genotypes)), 2))
    vars = [baseline_all,dF_F0_all]
    varnames = ['Baseline','dF_F0']

    for v,varname in enumerate(varnames):
        for i,__ in enumerate(combos):
            group1 = genotypes[combos[i][0]]
            group2 = genotypes[combos[i][1]]
            test = stats.mannwhitneyu(vars[v][combos[i][0]],vars[v][combos[i][1]])
            statistic = test.statistic
            pvalue = test.pvalue
            padj = pvalue*len(combos)

            result.append([varname,group1,group2,statistic,pvalue,padj])

    with open('{}/mannwhitneyu.txt'.format(output), 'w') as fp:
        for item in result:
            # write each item on a new line
            fp.write('%s\n' % ['group1','group2','statistic','pvalue','padjusted'])
            fp.write('%s\n' % item)
        print('Writing statistical test results...')
        
    baseline_mean = [np.mean(np.array(baseline_all[i])) for i in range(len(baseline_all))]
    baseline_sem = [stats.sem(np.array(baseline_all[i])) for i in range(len(baseline_all))]

    dF_F0_mean = [np.mean(np.array(dF_F0_all[i])) for i in range(len(dF_F0_all))]
    dF_F0_sem = [stats.sem(np.array(dF_F0_all[i])) for i in range(len(dF_F0_all))]


    fig,axs= plt.subplots(1,2, figsize=(3, 1.5), dpi=300)
    w = 0.8    # bar width
    x = np.arange(len(genotypes))

    axs[0].bar(x, height=baseline_mean, yerr=baseline_sem,capsize=3,width=w,color=colors)
    axs[1].bar(x, height=dF_F0_mean, yerr=dF_F0_sem,capsize=3,width=w,color=colors)

    for i in range(len(genotypes)):
    # distribute scatter randomly across whole width of bar
        axs[0].scatter(x[i] + np.random.random(len(baseline_all[i])) * w/2 - w/4, baseline_all[i], s = 3, c='black')
        axs[1].scatter(x[i] + np.random.random(len(dF_F0_all[i])) * w/2 - w/4, dF_F0_all[i], s = 3, c='black')


    axs[0].set_ylabel('Baseline F0')
    axs[1].set_ylabel('Mean dF/F0')
    #axs[1].set_ylim([0,3])

    for i in range(2):
        x_left, x_right = axs[i].get_xlim()
        y_low, y_high = axs[i].get_ylim()
        axs[i].set_aspect(abs((x_right-x_left)/(y_low-y_high))*2)
        axs[i].spines.right.set_visible(False)
        axs[i].spines.top.set_visible(False)

    plt.savefig('{}/barplot.svg'.format(output), format='svg')
    plt.savefig('{}/barplot.eps'.format(output), format='eps')

def analyze_excludeNoResponse(input,output,genotypes,colors,testType):
    import pickle
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    from itertools import combinations
    import scipy.stats as stats
    import scikit_posthocs as sp

    output = '{}/exclude_no_response'.format(output)
    if not os.path.exists(output):
        os.makedirs(output)

    data = []
    baseline_all = []
    dF_F0_all = []
    observations = []
    metrics = ['Raw ROI','dF','dF_F0']

    totalTime = 15

    for genotype in genotypes:
        raw,dF,dF_f0,mean_baseline,mean_dF_f0,max_dF_f0 = pickle.load(open('{}/{}/Values.pkl'.format(input,genotype), 'rb'))

        # empty list of boolean to store response
        #onResponse = []

        # first convert the reference dataframe to array
        #raw = np.array(raw)
        # total number of trials
        #n_larvae = len(raw)

        #for l in range(n_larvae):
            #if np.mean([raw[l][5:9]]) > (1.05*np.mean([raw[l][0:4]])) :
                #onResponse.append(True)
            #else:
                #onResponse.append(True)

        # drop values from trials with no responses or nan values
        #mean_baseline = [mean_baseline[l] for l in range(n_larvae) if (onResponse[l]) & (mean_baseline[l]!='nan')]
        #mean_dF_f0 = [mean_dF_f0[l] for l in range(n_larvae) if (onResponse[l]) & (mean_dF_f0[l]!='nan')]

        # compute response rate
        #observation = [sum(onResponse), n_larvae-sum(onResponse)]

        # to compute mean and sem of time series
        
        #vars = [raw,dF,dF_f0]

        # for each metric
        #for j,var in enumerate(vars):

            #mean = []
            #sem = []

            # subset var to include only responses
            #var = [var[l] for l in range(n_larvae) if onResponse[l]]
            # convert the variable to array
            #var = np.array(var)
            

            # then for each time point, compute the mean and sem over all larvae
            #for t in range(totalTime-1):
                #mean.append(np.mean([var[l][t] for l in range(len(var))]))
                #sem.append(stats.sem([var[l][t] for l in range(len(var))]))
            
            # append to the dataframe
            #data.append([mean,sem])

        baseline_all.append(mean_baseline)
        dF_F0_all.append(mean_dF_f0)
        #observations.append(observation)
    
    # convert the data to array    
    #data = np.array(data)
    # reshape the data
    #dims = data.shape
    #data = data.reshape(len(genotypes),len(metrics),dims[-2],dims[-1])

    # data is then a 3D list
    # 1st dim: genotype or conditions
    # 2nd dim: always of size 3, raw, dF, and dF_f0
    # 3rd dim: always of size 2, mean and sem

    # visualize
    # Say, "the default sans-serif font is COMIC SANS"
    plt.rcParams['font.sans-serif'] = "Arial"
    # Then, "ALWAYS use sans-serif fonts"
    plt.rcParams['font.family'] = "sans-serif"

    #fig,axs = plt.subplots(1,3, figsize=(11,8), dpi=300)
    #x = np.linspace(0,totalTime-1,totalTime-1)
    #for i,metric in enumerate(metrics):
        #for j,__ in enumerate(genotypes):
            #axs[i].plot(x,data[j][i][0], colors[j])
            #axs[i].fill_between(x,np.array(data[j][i][0])-np.array(data[j][i][1]),np.array(data[j][i][0])+np.array(data[j][i][1]), facecolor = colors[j], alpha = 0.3)
            #axs[i].set_title(metric)
            #axs[i].set_xlabel('Seconds (s)')

           # x_left, x_right = axs[i].get_xlim()
           # y_low, y_high = axs[i].get_ylim()       
           # axs[i].set_aspect(abs((x_right-x_left)/(y_low-y_high))*1)
           # axs[i].spines.right.set_visible(False)
           # axs[i].spines.top.set_visible(False)
    #plt.savefig('{}/lineplot.svg'.format(output), format='svg')
    #plt.savefig('{}/lineplot.pdf'.format(output), format='pdf')
    


    result = []
    combos = list(combinations(range(len(genotypes)), 2))
    vars = [baseline_all,dF_F0_all]
    varnames = ['Baseline','dF_F0']

    if testType == 'mannwhitneyu':
        # statistical test on mean baseline values and mean dF/F0 values
        for v,varname in enumerate(varnames):
            for i,__ in enumerate(combos):
                group1 = genotypes[combos[i][0]]
                group2 = genotypes[combos[i][1]]
                test = stats.mannwhitneyu(vars[v][combos[i][0]],vars[v][combos[i][1]])
                statistic = test.statistic
                pvalue = test.pvalue
                padj = pvalue*len(combos)

                result.append([varname,group1,group2,statistic,pvalue,padj])

        with open('{}/mannwhitneyu.txt'.format(output), 'w') as fp:
            for item in result:
                # write each item on a new line
                fp.write('%s\n' % ['group1','group2','statistic','pvalue','padjusted'])
                fp.write('%s\n' % item)
            print('Writing Mann-Whitney statistical test results...')
    
    if testType == 'kruskaldunn':
        pvals = []
        result = []

        with open('{}/kruskaldunn.txt'.format(output), 'w') as fp:
            print('Writing Kruskal-Wallis H-test and post-hoc Dunn test results...')
            fp.write('%s\n' % ['metric','genotypes','statistic','pvalue'])

            # statistical test on mean baseline values and mean dF/F0 values
            for v,varname in enumerate(varnames):
                # kruskal wallis test: non-parametric version of ANOVA to see global difference
                test = stats.kruskal(vars[v][0],vars[v][1],vars[v][2])
                statistic = test.statistic
                pvalue = test.pvalue
                result = [varname,genotypes,statistic,pvalue]
                # post-hoc Dunn test to see which are different
                pvals = sp.posthoc_dunn(vars[v])
                
                fp.write('%s\n' % result)
                fp.write('%s\n' % pvals)


    # statistical test on response rate
    #rate_result = []
    #for i,__ in enumerate(combos):
       # group1 = genotypes[combos[i][0]]
        #group2 = genotypes[combos[i][1]]

        #contigency_table = [
        #    observations[combos[i][0]],
        #    observations[combos[i][1]]
        #]

        # implement a fisher-exact test
        #from scipy.stats import fisher_exact
        #test = fisher_exact(contigency_table)

        #statistic = test[0]
        #pvalue = test[1]
        #padj = pvalue*len(combos)
        #rate_result.append([group1,group2,contigency_table,statistic,pvalue,padj,])

    #with open('{}/fisher_exact.txt'.format(output), 'w') as fp:
        #for item in rate_result:
            # write each item on a new line
            #fp.write('%s\n' % ['group1','group2','contigency','statistic','pvalue','padjusted'])
            #fp.write('%s\n' % item)
        #print('Writing Fisher exact statistical test results...')

        
    baseline_mean = [np.mean(np.array(baseline_all[i])) for i in range(len(baseline_all))]
    baseline_sem = [stats.sem(np.array(baseline_all[i])) for i in range(len(baseline_all))]

    dF_F0_mean = [np.mean(np.array(dF_F0_all[i])) for i in range(len(dF_F0_all))]
    dF_F0_sem = [stats.sem(np.array(dF_F0_all[i])) for i in range(len(dF_F0_all))]

########### bar plots ###############
    fig,axs= plt.subplots(1,2, figsize=(4,1.5), dpi=300)
    w = 0.8    # bar width
    x = np.arange(len(genotypes))

    axs[0].bar(x, height=baseline_mean, yerr=baseline_sem,capsize=3,width=w,color=colors)
    axs[1].bar(x, height=dF_F0_mean, yerr=dF_F0_sem,capsize=3,width=w,color=colors)

    #for i in range(len(genotypes)):
    # distribute scatter randomly across whole width of bar
    #    axs[0].scatter(x[i] + np.random.random(len(baseline_all[i])) * w/2 - w/4, baseline_all[i], s = 3, c='black')
    #    axs[1].scatter(x[i] + np.random.random(len(dF_F0_all[i])) * w/2 - w/4, dF_F0_all[i], s = 3, c='black')


    axs[0].set_ylabel('Baseline F0')
    axs[1].set_ylabel('Mean dF/F0')
    #axs[1].set_ylim([0,3])

    for i in range(2):
        x_left, x_right = axs[i].get_xlim()
        y_low, y_high = axs[i].get_ylim()
        axs[i].set_aspect(abs((x_right-x_left)/(y_low-y_high))*2)
        axs[i].spines.right.set_visible(False)
        axs[i].spines.top.set_visible(False)
        axs[i].set_xticks([])
        axs[i].xaxis.set_tick_params(labelbottom=False)

    plt.savefig('{}/barplot.svg'.format(output), format='svg')
    plt.savefig('{}/barplot.eps'.format(output), format='eps')

########### violin plots ###############
    fig,axs= plt.subplots(1,2, figsize=(4,1.5), dpi=300)
    w = 0.8    # bar width
    x = np.arange(len(genotypes))

    bplot1 = axs[0].boxplot(baseline_all, widths = 0.7*w, patch_artist = True, medianprops = dict(color = "black"), flierprops={'marker': 'o', 'markersize': 2})
    bplot2 = axs[1].boxplot(dF_F0_all, widths = 0.7*w, patch_artist = True, medianprops = dict(color = "black"),flierprops={'marker': 'o', 'markersize': 2})


    for bplot in (bplot1, bplot2):
        for patch, color in zip(bplot['boxes'], colors):
            patch.set_facecolor(color)

    axs[0].set_ylabel('Baseline F0')
    axs[1].set_ylabel('Mean dF/F0')

    for i in range(2):
        x_left, x_right = axs[i].get_xlim()
        y_low, y_high = axs[i].get_ylim()
        axs[i].set_aspect(abs((x_right-x_left)/(y_low-y_high))*2)
        axs[i].spines.right.set_visible(False)
        axs[i].spines.top.set_visible(False)
        axs[i].set_xticks([])
        axs[i].xaxis.set_tick_params(labelbottom=False)

    plt.savefig('{}/boxplot.svg'.format(output), format='svg')
    plt.savefig('{}/boxplot.eps'.format(output), format='eps')

    # plot response rate
    #fig,ax= plt.subplots(1,1, figsize=(3, 1.5), dpi=300)
    #ax.bar(x,[observations[j][0]/sum(observations[j]) for j in range(len(genotypes))],color=colors)
  
    #ax.set_ylabel('Response Probability')
    #ax.set_ylim([0,1])

    #x_left, x_right = ax.get_xlim()
    #y_low, y_high = ax.get_ylim()       
    #ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*1.5)
    #ax.spines.right.set_visible(False)
    #ax.spines.top.set_visible(False)
    #ax.set_xticks([])
    #ax.xaxis.set_tick_params(labelbottom=False)

    #plt.savefig('{}/responseRate.svg'.format(output), format='svg')
    #plt.savefig('{}/responseRate.eps'.format(output), format='eps')

def analyze_cdf(input,output,genotypes,colors):
    import pickle
    import numpy as np
    import os
    import matplotlib.pyplot as plt
    from itertools import combinations
    import scipy.stats as stats

    output = '{}/cdf'.format(output)
    if not os.path.exists(output):
        os.makedirs(output)

    data = []
    baseline_all = []
    mean_dF_F0_all = []
    max_dF_F0_all = []
    metrics = ['Raw ROI','dF','dF_F0']

    totalTime = 15

    for genotype in genotypes:
        raw,dF,dF_f0,mean_baseline,mean_dF_f0, max_dF_f0 = pickle.load(open('{}/{}/Values.pkl'.format(input,genotype), 'rb'))
        
        vars = [raw,dF,dF_f0]

        # for each metric
        for j,var in enumerate(vars):

            mean = []
            sem = []
            # convert the variable to array
            var = np.array(var)

            # then for each time point, compute the mean and sem over all larvae
            for t in range(totalTime-1):
                mean.append(np.mean([var[l][t] for l in range(len(var))]))
                sem.append(stats.sem([var[l][t] for l in range(len(var))]))
            
            # append to the dataframe
            data.append([mean,sem])

        
        baseline_all.append(mean_baseline)
        mean_dF_F0_all.append(mean_dF_f0)
        max_dF_F0_all.append(max_dF_f0)

        
    # convert the data to array    
    data = np.array(data)
    # reshape the data
    dims = data.shape
    data = data.reshape(len(genotypes),len(metrics),dims[-2],dims[-1])

    # data is then a 3D list
    # 1st dim: genotype or conditions
    # 2nd dim: always of size 3, raw, dF, and dF_f0
    # 3rd dim: always of size 2, mean and sem

    # visualize
    # Say, "the default sans-serif font is COMIC SANS"
    plt.rcParams['font.sans-serif'] = "Arial"
    # Then, "ALWAYS use sans-serif fonts"
    plt.rcParams['font.family'] = "sans-serif"

    fig,axs = plt.subplots(1,3, figsize=(11,8), dpi=300)
    x = np.linspace(0,totalTime-1,totalTime-1)
    for i,metric in enumerate(metrics):
        for j,__ in enumerate(genotypes):
            axs[i].plot(x,data[j][i][0], colors[j])
            axs[i].fill_between(x,np.array(data[j][i][0])-np.array(data[j][i][1]),np.array(data[j][i][0])+np.array(data[j][i][1]), facecolor = colors[j], alpha = 0.3)
            axs[i].set_title(metric)
            axs[i].set_xlabel('Seconds (s)')

            x_left, x_right = axs[i].get_xlim()
            y_low, y_high = axs[i].get_ylim()       
            axs[i].set_aspect(abs((x_right-x_left)/(y_low-y_high))*1)
            axs[i].spines.right.set_visible(False)
            axs[i].spines.top.set_visible(False)
    plt.savefig('{}/lineplot.svg'.format(output), format='svg')
    plt.savefig('{}/lineplot.pdf'.format(output), format='pdf')


    result = []
    combos = list(combinations(range(len(genotypes)), 2))
    vars = [mean_dF_F0_all,max_dF_F0_all,baseline_all]
    varnames = ['Mean dF/F0','Maximum dF/F0','Mean Baseline Fluorescence']

    for v,varname in enumerate(varnames):
        for i,__ in enumerate(combos):
            group1 = genotypes[combos[i][0]]
            group2 = genotypes[combos[i][1]]
            test = stats.kstest(vars[v][combos[i][0]],vars[v][combos[i][1]])
            statistic = test.statistic
            pvalue = test.pvalue
            # adjust using Bonferroni's correction
            padj = pvalue*len(combos)

            result.append([varname,group1,group2,statistic,pvalue,padj])

    with open('{}/Kolmogorov-Smirnov.txt'.format(output), 'w') as fp:
        for item in result:
            # write each item on a new line
            fp.write('%s\n' % ['group1','group2','statistic','pvalue','padjusted'])
            fp.write('%s\n' % item)
        print('Writing statistical test results...')
    
    # visualize the cdf
    plt.rcParams['font.sans-serif'] = "Arial"
    # Then, "ALWAYS use sans-serif fonts"
    plt.rcParams['font.family'] = "sans-serif"

    fig,axs = plt.subplots(1,len(vars), figsize=(6,3), dpi=300)
    
    for i,__ in enumerate(vars):
        for j,__ in enumerate(genotypes):
            count, bins_count = np.histogram(vars[i][j], bins=50)
            # finding the PDF of the histogram using count values
            pdf = count/sum(count)
            # using numpy np.cumsum to calculate the CDF
            # # We can also find using the PDF values by looping and adding
            cdf = np.cumsum(pdf)
            # plot the cdf
            axs[i].step(bins_count[1:], cdf, color=colors[j])
            axs[i].set_xlabel(varnames[i])
            axs[i].set_ylabel('CDF')
            x_left, x_right = axs[i].get_xlim()
            y_low, y_high = axs[i].get_ylim()       
            axs[i].set_aspect(abs((x_right-x_left)/(y_low-y_high))*1)
            axs[i].spines.right.set_visible(False)
            axs[i].spines.top.set_visible(False)
    plt.tight_layout()
    plt.savefig('{}/cdf.eps'.format(output), format='eps')

def correlate(root_folder,paths,output):
    import glob
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    totalTime = 15

    # to store info
    baselines = []
    mean_dFs = []
    genotypes = []
    trials = []

    # for each genotype given
    for path in paths:

        # find the corresponding folder that contains raw data for that specific trial
        files = glob.glob("{}/{}/**/Values.csv".format(root_folder,path), recursive = True)

        # for each trial
        for file in files:
            # append the trial number extracted from file name
            trials.append(file.split('/')[-2].split('_')[-2])
            
            with open(file) as f:
                # open Values.csv file, then process
                lines = (line for line in f if not line.startswith('#'))
                data = np.loadtxt(lines, delimiter=',', skiprows=1)
                data = np.transpose(data)   
                data = data[-1]

                if len(data) > totalTime:
                    data = data[0:totalTime-1]

                # mean baseline value is the average of the first 5 elements
                baseline = np.mean(data[0:4]) 
                baselines.append(baseline)
                # mean dF is the average of the next 5 elements
                mean_dF = (data-baseline)/baseline
                mean_dFs.append(np.mean(mean_dF[5:9]))
                # note the genotype of trial, the number of trial, and the strength of stimuli
                genotypes.append(path)

    # zip together as pd dataframe                
    df = pd.DataFrame(list(zip(genotypes, trials, baselines, mean_dFs)),
        columns =['genotype', 'trial','baseline','mean_dF_F0'])
    # convert the trial factor from string to integer
    df['trial'] = df['trial'].astype(int)
    df = df[df['trial']<=8]

    import seaborn as sns
    # global figure setting
    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.family'] = "sans-serif"

    for genotype in df['genotype'].unique():
        # print the correlation
        print(df[df['genotype']==genotype].corr())

        # visualize using violin plot
        fig,axs = plt.subplots(1,2, figsize=(8,3), dpi=300)

        sns.violinplot(data=df[df['genotype']==genotype], x="trial", y="baseline", ax = axs[0])
        # plot the median value of first trial to compare to the rest
        axs[0].axhline(y=df[(df['trial']==1) & (df['genotype']==genotype)]['baseline'].median(), color='r', linestyle='-')
        axs[0].set_title(genotype)

        axs[1] = sns.violinplot(data=df[df['genotype']==genotype], x="trial", y="mean_dF_F0",ax = axs[1])
        # plot the median value of first trial to compare to the rest
        axs[1].axhline(y=df[(df['trial']==1) & (df['genotype']==genotype)]['mean_dF_F0'].median(), color='r', linestyle='-')
        axs[1].set_title(genotype)
        genotype = str(genotype).replace('/','_')
        plt.savefig('{}/{}_violinplot.eps'.format(output,genotype), format='eps')

def heatmap(input,output,genotypes):
    import pickle
    import numpy as np
    import os
    import matplotlib.pyplot as plt

    totalTime = 15

    output = '{}/heatmaps'.format(output)
    if not os.path.exists(output):
        os.makedirs(output)

    metrics = ['Raw ROI','dF','dF_F0']

    maximum = []
    minimum = []

    for genotype in genotypes:
        raw,dF,dF_f0,mean_baseline,mean_dF_f0,max_dF_f0 = pickle.load(open('{}/{}/Values.pkl'.format(input,genotype), 'rb'))
        # compute the range of values wanted
        maximum.append([max_value(raw),max_value(dF),max_value(dF_f0)])
        minimum.append([min_value(raw),min_value(dF),min_value(dF_f0)])
    
    # transpose the matrix so that the row dimension is of the same type of metrics over all genotypes
    maximum = np.array(maximum).T
    minimum = np.array(minimum).T
    # find the maximum value of each metric across all genotypes and store in 1D list
    maximum = [max(item) for item in maximum]
    minimum = [min(item) for item in minimum]

    # plot the heatmaps
    for genotype in genotypes:

        raw,dF,dF_f0,mean_baseline,mean_dF_f0,max_dF_f0 = pickle.load(open('{}/{}/Values.pkl'.format(input,genotype), 'rb'))

        plt.rcParams['font.sans-serif'] = "Arial"
        plt.rcParams['font.family'] = "sans-serif"

        fig,axs = plt.subplots(3,len([raw,dF,dF_f0]), figsize = (9,9), dpi = 300)

        for i,item in enumerate([raw,dF,dF_f0]):

            item = np.vstack([row[0:totalTime-1] for row in item])

            means = []
            means_off = []
            for j in range(len(item)):
                means.append(np.mean(item[j,5:9]))
                means_off.append(np.mean(item[j,10:14]))
            sort_index = np.argsort(means)
            item_sorted_stim = np.vstack([item[j] for j in sort_index])
            sort_index = np.argsort(means_off)
            item_sorted_off = np.vstack([item[j] for j in sort_index])

            im = axs[0,i].imshow(item, cmap='jet', interpolation='nearest', vmin = minimum[i],vmax = maximum[i])
            plt.colorbar(im, ax=axs[0,i])
            axs[0,i].set_title(metrics[i])
            axs[0,i].set_xlabel('Seconds (s)')
            x_left, x_right = axs[0,i].get_xlim()
            y_low, y_high = axs[0,i].get_ylim()       
            axs[0,i].set_aspect(abs((x_right-x_left)/(y_low-y_high))*1)

            im = axs[1,i].imshow(item_sorted_stim, cmap='jet', interpolation='nearest',vmin = minimum[i],vmax = maximum[i])
            plt.colorbar(im, ax=axs[1,i])
            axs[1,i].set_title(metrics[i])
            axs[1,i].set_xlabel('Seconds (s)')
            x_left, x_right = axs[1,i].get_xlim()
            y_low, y_high = axs[1,i].get_ylim()       
            axs[1,i].set_aspect(abs((x_right-x_left)/(y_low-y_high))*1)

            im = axs[2,i].imshow(item_sorted_off, cmap='jet', interpolation='nearest',vmin = minimum[i],vmax = maximum[i])
            plt.colorbar(im, ax=axs[2,i])
            axs[2,i].set_title(metrics[i])
            axs[2,i].set_xlabel('Seconds (s)')
            x_left, x_right = axs[2,i].get_xlim()
            y_low, y_high = axs[2,i].get_ylim()       
            axs[2,i].set_aspect(abs((x_right-x_left)/(y_low-y_high))*1)

            plt.savefig('{}/{}.eps'.format(output,str(genotype).split('/')[0]), format='eps')

    
def average_each_larva(root_folder,paths,output):
    import glob
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    totalTime = 15

    # to store info
    baselines = []
    mean_dFs = []
    genotypes = []
    larvae = []
    dates = []

    # for each genotype given
    for path in paths:

        # find the corresponding folder that contains raw data for that specific trial
        files = glob.glob("{}/{}/**/Values.csv".format(root_folder,path), recursive = True)

        # for each trial
        for file in files:
            # append the larva id and date of experiment extracted from file name
            larvae.append(file.split('/')[-2].split('_')[-3].split('-')[0])
            dates.append(file.split('/')[-2].split('_')[1])
            
            with open(file) as f:
                # open Values.csv file, then process
                lines = (line for line in f if not line.startswith('#'))
                data = np.loadtxt(lines, delimiter=',', skiprows=1)
                data = np.transpose(data)   
                data = data[-1]

                if len(data) > totalTime:
                    data = data[0:totalTime-1]

                # mean baseline value is the average of the first 5 elements
                baseline = np.mean(data[0:4]) 
                baselines.append(baseline)
                # mean dF is the average of the next 5 elements
                mean_dF = (data-baseline)/baseline
                mean_dFs.append(np.mean(mean_dF[5:9]))
                # note the genotype of trial, the number of trial, and the strength of stimuli
                genotypes.append(path)

    # zip together as pd dataframe                
    df = pd.DataFrame(list(zip(genotypes, larvae, dates, baselines, mean_dFs)),
        columns =['genotype', 'larva_number', 'date', 'baseline','mean_dF_F0'])
    # convert the trial factor from string to integer
    df['larva_number'] = df['larva_number'].astype(int)

    import seaborn as sns
    # global figure setting
    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.family'] = "sans-serif"

    for genotype in df['genotype'].unique():
        # print the correlation
        print(df[df['genotype']==genotype].corr())

        # visualize using violin plot
        fig,axs = plt.subplots(1,2, figsize=(8,3), dpi=300)

        axs[0] = sns.violinplot(data=df[df['genotype']==genotype], x="date", y="baseline", hue="larva_number", ax = axs[0],legend=False)
        # plot the median value of first trial to compare to the rest
        axs[0].set_title(genotype)
        axs[0].legend_.remove()

        axs[1] = sns.violinplot(data=df[df['genotype']==genotype], x="date", y="mean_dF_F0", hue = "larva_number", ax = axs[1],legend=False)
        # plot the median value of first trial to compare to the rest
        axs[1].set_title(genotype)
        axs[1].legend_.remove()

        genotype = str(genotype).replace('/','_')
        plt.savefig('{}/{}_averageLarva_violinplot.eps'.format(output,genotype), format='eps')

def ratio_ONvsOFF(input,output,paths):
    import glob
    import os
    import numpy as np
    import pandas as pd
    import pickle

    totalTime = 15

    responseFref = []
    responseProb = []
    genotypes = []

    for path in paths:
        #for each path
        files = glob.glob("{}/{}/**/Values.csv".format(input,path), recursive = True)

        ON = 0
        OFF = 0
        ON_OFF = 0
        no_response = 0

        for file in files:

            with open(file) as f:
                lines = (line for line in f if not line.startswith('#'))
                data = np.loadtxt(lines, delimiter=',', skiprows=1)
                data = np.transpose(data)   
                data = data[-1]

                # make sure data is of right length
                if len(data) > totalTime:
                    data = data[0:totalTime-1]

                # ON response = dF highest during stimulation
                if np.mean([data[5:9]]) > (1.05*np.mean([data[0:4]])) :
                    ON += 1 
                # OFF response = dF highest right after stimulation, no response during stimulation
                elif (np.mean([data[-5:]]) > (1.05*np.mean([data[0:4]]))) and (np.mean([data[5:9]]) <= (1.05*np.mean([data[0:4]]))):
                    OFF += 1
                # ON response and even higher OFF response
                elif (np.mean([data[5:9]]) > (1.05*np.mean([data[0:4]]))) and (np.mean([data[-5:]]) >= (1.05*np.mean([data[5:9]]))):
                    ON_OFF += 1
                else:
                    no_response += 1

        my_array = np.array([ON,OFF,ON_OFF,no_response])
        responseFref.append(my_array)
        responseProb.append(my_array/np.sum(my_array))
        genotypes.append(path) 
    
    output = '{}/ratio_ONvsOFF'.format(output)
    if not os.path.exists(output):
        os.makedirs(output)
    
    with open('{}/ratios.txt'.format(output), 'w') as fp:
        fp.write('%s\n' % genotypes)
        fp.write('%s\n' % responseFref)
        fp.write('%s\n' % responseProb)
    
    # transpose the matrix
    responseProbs = np.array(responseProb).T
    responseTypes = ['ON','OFF','ON_OFF','no_response']
    bars = {responseType:responseProbs[i] for i, responseType in enumerate(responseTypes)}

    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif'] = "Arial"
    # Then, "ALWAYS use sans-serif fonts"
    plt.rcParams['font.family'] = "sans-serif"

    fig, ax = plt.subplots(1,1, figsize=(3, 1.5), dpi=300)
    bottom = np.zeros(3)

    colors = ['black','#777777','#999999','#cccccc']
    for boolean, bar in bars.items():
        p = ax.bar(genotypes, bar, label = boolean, color = colors[list(bars.keys()).index(boolean)], bottom=bottom)
        bottom += bar

        ax.set_ylabel('Probability')
        #ax.set_ylim([0,1])
        x_left, x_right = ax.get_xlim()
        y_low, y_high = ax.get_ylim()       
        ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*1.5)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.set_xticks([])
        ax.xaxis.set_tick_params(labelbottom=False)
        plt.legend(loc=(1.04, 0))

    plt.savefig('{}/responseRatios.svg'.format(output), format='svg')
    plt.savefig('{}/responseRatios.eps'.format(output), format='eps')
