#!/usr/bin/python

from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

#   subroutines to analyze Pol-II ChIP experiments

def positives(df): # get positive indiced from df from limma
    return df[(df.logFC>np.log(1.5)) & (df['adj.P.Val']<=0.05)].index
    
def negatives(df): # get negative indices from df from limma
    return df[(df.logFC<-np.log(1.5)) & (df['adj.P.Val']<=0.05)].index

def volcano(df):
    # plot a volcano plot
    p = -np.log10(df['adj.P.Val'])
    plt.scatter(df.logFC, p, s=3, alpha=0.4, c='k')
    maxP = np.max(p)
    plt.plot([df.logFC.min(), df.logFC.max()], [-np.log10(0.05),-np.log10(0.05)], ls="--", c='r')
    plt.plot([-np.log(1.5), -np.log(1.5)], [0,maxP], ls="--", c='r')
    plt.plot([np.log(1.5), np.log(1.5)], [0,maxP], ls="--", c='r')
    #plt.show()
    #plt.close()
    
    
def doPCA(df, plot=True):  ## quizas escribir *args or **kargs for {scatter: [], text:[]}
    '''
        INPUT: df of counts. Rows are genes and columns are experiments (e.g. c1,c2,c3,t1,t2,t3...)
        OUTPUT: explained variance between samples in the 2 first principal components
    '''
    # first order the df for text to be correctly alligned
    ctrols = [i for i in df.columns if i[0]=='c']
    tests = [i for i in df.columns if i[0]=='t' ]
    df = df[ctrols+tests]
    pca1 = PCA(n_components=2).fit(df.values.T)
    pca = pca1.transform(df.values.T)

    if plot:        
        colors = ['r']*len(ctrols) + ['b']*len(tests)
        plt.scatter(pca[:,0], pca[:,1], c=colors)
        text = ctrols+tests
        [plt.text(pca[n,0]*1.1, pca[n,1]*1.1, text[n]) for n in range(len(df.columns))]
        plt.grid()

    return pca1.explained_variance_ratio_



class limado:
    def __init__(self, filename):
        self.filename = filename
        self.df = pd.read_csv(filename, index_col=0)
        self.dist2tss = self.df.loc[-50:250].idxmax()
        self.norm_counts = self.counts()/self.spike()
        self.profile = self.df.mean(axis=1)
        
    def counts(self):
        # counts of reads where the pol2 is transcribing at the beginning of each gene
        v, k = [], []
        for i in self.dist2tss.index:
            m = self.dist2tss.loc[i]
            v.append(self.df.loc[m-50:m+50, i].sum())
            k.append(i)
        return pd.DataFrame(v, index=k) 

    def spike(self):
        # number of reads of pombe to normalize the counts
        pombe = self.filename[:-8] + 'pombe.sam'
        return sum([1 for i in open(self.filename) if i[0]!='@'])

    
''' Directories where files obtained by
    1- Aligning to genome with Bowtie2
    2- annotating and counting with metGen.py
    are stored.
'''    
spt8_folder = ['1_Sample_161_Spt8_deg_YPD_DMSO_A/',
               '3_Sample_163_Spt8_deg_YPD_DMSO_B/',
               '2_Sample_162_Spt8_deg_YPD_3IAA_A/',
               '4_Sample_164_Spt8_deg_YPD_3IAA_B/']

spt7_3_folder = ['10_Sample_252A_Spt7_3_YPD_DMSO_B/',
                 '11_Sample_252B_Spt7_3_YPD_DMSO_B/',
                 '12_Sample_252C_Spt7_3_YPD_DMSO_B/',
                 '1_Sample_249A_Spt7_3_YPD_3IAA_A/',
                 '2_Sample_249B_Spt7_3_YPD_3IAA_A/',
                 '3_Sample_249C_Spt7_3_YPD_3IAA_A/',
                 '4_Sample_250A_Spt7_3_YPD_DMSO_A/',
                 '5_Sample_250B_Spt7_3_YPD_DMSO_A/',
                 '6_Sample_250C_Spt7_3_YPD_DMSO_A/',
                 '7_Sample_251A_Spt7_3_YPD_3IAA_B/',
                 '8_Sample_251B_Spt7_3_YPD_3IAA_B/',
                 '9_Sample_251C_Spt7_3_YPD_3IAA_B/']

spt20_ypd_folder = ['10_Sample_144_Spt20_deg_YPD_3IAA_B/',
                     '15_Sample_221_Spt20_deg_YPD_DMSO_A/',
                     '16_Sample_222_Spt20_deg_YPD_3IAA_A/',
                     '17_Sample_223_Spt20_deg_YPD_DMSO_B/',
                     '18_Sample_224_Spt20_deg_YPD_3IAA_B/',
                     '8_Sample_142_Spt20_deg_YPD_3IAA_A/',
                     '9_Sample_143_Spt20_deg_YPD_DMSO_B/']

spt20_gc_folder = ['11_Sample_145_Spt20_deg_GC_DMSO_A/',
                   '12_Sample_146_Spt20_deg_GC_3IAA_A/',
                   '13_Sample_147_Spt20_deg_GC_DMSO_B/',
                   '14_Sample_148_Spt20_deg_GC_3IAA_B/',]

spt8_3_folder = ['10_Sample_254C_Spt8_3_YPD_DMSO_A/',
                  '11_Sample_255A_Spt8_3_YPD_3IAA_B/',
                  '12_Sample_255B_Spt8_3_YPD_3IAA_B/',
                  '13_Sample_255C_Spt8_3_YPD_3IAA_B/',
                  '14_Sample_256A_Spt8_3_YPD_DMSO_B/',
                  '15_Sample_256B_Spt8_3_YPD_DMSO_B/',
                  '16_Sample_256C_Spt8_3_YPD_DMSO_B/',
                  '1_Sample_217_Spt3_Spt8_deg_YPD_DMSO_A/',
                  '2_Sample_218_Spt3_Spt8_deg_YPD_3IAA_A/',
                  '3_Sample_219_Spt3_Spt8_deg_YPD_DMSO_B/',
                  '4_Sample_220_Spt3_Spt8_deg_YPD_3IAA_B/',
                  '5_Sample_253A_Spt8_3_YPD_3IAA_A/',
                  '6_Sample_253B_Spt8_3_YPD_3IAA_A/',
                  '7_Sample_253C_Spt8_3_YPD_3IAA_A/',
                  '8_Sample_254A_Spt8_3_YPD_DMSO_A/',
                  '9_Sample_254B_Spt8_3_YPD_DMSO_A/']

spt20_3_folder = ['10_Sample_248A_Spt20_3_YPD_DMSO_B/',
                   '11_Sample_248B_Spt20_3_YPD_DMSO_B/',
                   '12_Sample_248C_Spt20_3_YPD_DMSO_B/',
                   '1_Sample_245A_Spt20_3_YPD_3IAA_A/',
                   '2_Sample_245B_Spt20_3_YPD_3IAA_A/',                                                                                                                            '3_Sample_245C_Spt20_3_YPD_3IAA_A/',
                   '4_Sample_246A_Spt20_3_YPD_DMSO_A/',
                   '5_Sample_246B_Spt20_3_YPD_DMSO_A/',
                   '6_Sample_246C_Spt20_3_YPD_DMSO_A/',
                   '7_Sample_247A_Spt20_3_YPD_3IAA_B/',
                   '8_Sample_247B_Spt20_3_YPD_3IAA_B/',
                   '9_Sample_247C_Spt20_3_YPD_3IAA_B/']
    
spt20_7_folder = ['10_Sample_244A_Spt20_7_YPD_DMSO_B/',                                                                   
                  '11_Sample_244B_Spt20_7_YPD_DMSO_B/',
                  '12_Sample_244C_Spt20_7_YPD_DMSO_B/',
                  '1_Sample_241A_Spt20_7_YPD_3IAA_A/',
                  '2_Sample_241B_Spt20_7_YPD_3IAA_A/',
                  '3_Sample_241C_Spt20_7_YPD_3IAA_A/',
                  '4_Sample_242A_Spt20_7_YPD_DMSO_A/',
                  '5_Sample_242B_Spt20_7_YPD_DMSO_A/',
                  '6_Sample_242C_Spt20_7_YPD_DMSO_A/',
                  '7_Sample_243A_Spt20_7_YPD_3IAA_B/',
                  '8_Sample_243B_Spt20_7_YPD_3IAA_B/',
                  '9_Sample_243C_Spt20_7_YPD_3IAA_B/']


    
#def spike(filename):
#    # I might end up using normalization by spike in ONLY...
#    pombe = filename[:-8] + 'pombe.sam'
#    return sum([1 for i in open(filename) if i[0]!='@'])
#
#def limado(filename):
#    '''
#        This is the data that I will load on limma for the empirical bayes.
#        input = filename
#        output = {limma: values of peak+-50bp
#                  profile: profile of TSS+-1000bp
#                  M: distance to TSS of every gene}
#    '''
#    df = pd.read_csv(filename, index_col=0)
#    M = df.loc[-50:250].idxmax()
#    v, k = [], []
#    for i in M.index:
#        m = M.loc[i]
#        v.append(df.loc[m-50:m+50, i].sum())
#        k.append(i)
#    # finally also make the profile to have it handy
#    profile = df.mean(axis=1)
#    limma = pd.DataFrame(v, index=k) / spike(filename)
#    return limma, profile, M