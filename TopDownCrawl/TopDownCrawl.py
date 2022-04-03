#! /usr/bin/env python

import sys
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import logomaker as lm
from collections import defaultdict
from fonts.ttf import Roboto

from matplotlib import font_manager
font_manager.fontManager.addfont(Roboto)

N = ["A", "C", "G", "T"]
NN = ['AA', 'CA', 'GA', 'TA', 'AC', 'CC', 'GC', 'TC', 'AG', 'CG', 'GG', 'TG', 'AT', 'CT', 'GT', 'TT']
alts = {"A":["C", "G", "T"],"C":["A", "G", "T"],"G":["A", "C", "T"],"T":["A", "C", "G"]}

comp = {"A":"T", "C":"G", "G":"C", "T":"A", "N":"N", "_":"_"}
def rc(seqs):
    return np.array(["".join([comp[x] for x in seq][::-1]) for seq in seqs])

def fill(i, shift, crawl):
    crawl["Shift"][i] = shift
    crawl["Filled"][i] = True
    return crawl

def crawlSNPs(i, crawl, lookup):
    seq = crawl["Seqs"][i]
    for j,c in enumerate(seq):
        for d in alts[c]:
            temp = seq[0:j] + d + seq[j+1:]
            desti = lookup.pop(temp, -1)
            rci = lookup.pop(rc([temp])[0], -1)
            if(not desti == -1):
                crawl["Checked"][rci] = True
                crawl = fill(desti, crawl["Shift"][i], crawl)
    return crawl

def crawlLeft(i, crawl, lookup):
    for c in N:
        temp = c + crawl["Seqs"][i][:-1]
        desti = lookup.pop(temp, -1)
        rci = lookup.pop(rc([temp])[0], -1)
        if(not desti == -1):
            crawl["Checked"][rci] = True
            crawl = fill(desti, crawl["Shift"][i]-1, crawl)
            # crawl = crawlSNPs(desti, crawl, lookup, temp)
    for c in NN:
        temp = c + crawl["Seqs"][i][:-2]
        desti = lookup.pop(temp, -1)
        rci = lookup.pop(rc([temp])[0], -1)
        if(not desti == -1):
            crawl["Checked"][rci] = True
            crawl = fill(desti, crawl["Shift"][i]-2, crawl)
    return crawl
    
def crawlRight(i, crawl, lookup):
    for c in N:
        temp = crawl["Seqs"][i][1:] + c
        desti = lookup.pop(temp, -1)
        rci = lookup.pop(rc([temp])[0], -1)
        if(not desti == -1):
            crawl["Checked"][rci] = True
            crawl = fill(desti, crawl["Shift"][i]+1, crawl)
            # crawl = crawlSNPs(desti, crawl, lookup, temp)
    for c in NN:
        temp = crawl["Seqs"][i][2:] + c
        desti = lookup.pop(temp, -1)
        rci = lookup.pop(rc([temp])[0], -1)
        if(not desti == -1):
            crawl["Checked"][rci] = True
            crawl = fill(desti, crawl["Shift"][i]+2, crawl)
    return crawl

def plotPWM(imagename, seqs, weights):
    weights = np.array(weights).reshape(-1,1)
    seqs = np.array([list(seq) for seq in seqs])
    X = [(seqs == c).astype(float) * weights for c in N]
    X = [np.sum(x, axis=0) for x in X]
    X = np.array(X).transpose()
    X = pd.DataFrame(X, columns=N)
    X = lm.transform_matrix(X, from_type='counts', to_type='information')
    logo = lm.Logo(X, font_name="Roboto", 
        color_scheme={
            'A': (16/255, 150/255, 72/255),
            'C':(37/255, 92/255, 153/255),
            'G':(247/255, 179/255, 43/255),
            'T':(214/255, 40/255, 57/255)})
    logo.ax.set_xlabel('Position',fontsize=14)
    logo.ax.set_ylabel("Bits", fontsize=14)
    plt.savefig(imagename, bbox_inches="tight", dpi=600)
    plt.close()
    
def TDC(filename):
    df = None
    if(filename.endswith(".xlsx") or filename.endswith(".xls")):
        df = pd.read_excel(filename, usecols=(0, 1))
    elif(filename.endswith(".csv")):
        df = pd.read_csv(filename, sep=",", usecols=(0, 1))
    else:
        df = pd.read_csv(filename, delim_whitespace=True, usecols=(0, 1))
    
    basename = ".".join(filename.split(".")[:-1])
    scorename = df.columns[1]
    df_rc = df.copy()
    df_rc.iloc[:,0] = rc(df.iloc[:,0])
    df = pd.concat((df, df_rc)).iloc[:,0:2]
    df = df.groupby(df.columns[0]).mean().sort_values(scorename, ascending=False).reset_index()
    df = df.rename(columns={scorename:scorename + "_avg"})
    scorename = scorename + "_avg"

    seqs = df.iloc[:,0].values
    scores = df.iloc[:,1].values

    crawl = {"Seqs":seqs, scorename:scores, "Shift":np.zeros(len(seqs), dtype=int),
        "Filled":np.full(len(seqs), False), "Checked":np.full(len(seqs), False)}

    lookup = defaultdict(lambda: -1, zip(crawl["Seqs"], range(len(crawl["Seqs"]))))
    crawl = fill(0,0,crawl)
    # print("Aligning", len(lookup), "sequences . . .")
    lookup.pop(crawl["Seqs"][0])
    if(not crawl["Seqs"][0] == rc([crawl["Seqs"][0]])[0]):
        crawl["Checked"][lookup.pop(rc([crawl["Seqs"][0]])[0])] = True
    nextMax = np.argmax(~crawl["Checked"] * crawl["Filled"])
    while(crawl["Checked"][nextMax] == False and len(lookup) > 0):
        crawl["Checked"][nextMax] = True
        crawl = crawlSNPs(nextMax, crawl, lookup)
        crawl = crawlLeft(nextMax, crawl, lookup)
        crawl = crawlRight(nextMax, crawl, lookup)
        nextMax = np.argmax(~crawl["Checked"] * crawl["Filled"])
    # print("Alignment complete")
    print("Unable to align",len(lookup),"sequences")


    crawl = pd.DataFrame(crawl)
    crawl = crawl[crawl["Filled"] == True]

    crawl = crawl[["Seqs", scorename, "Shift"]]

    minShift = -np.min(crawl["Shift"])
    maxShift = np.max(crawl["Shift"])
    for i, item in crawl.iterrows():
        shift = item[2]
        crawl.at[i, "Seqs"] = "_" * (minShift+shift) + item[0] + "_" * (maxShift-shift)

    crawl.to_csv(basename + "_aligned.tsv", index=False, header=crawl.columns, sep="\t")

    ltrim = max(0, min(2, minShift-2))
    rtrim = max(0, min(2, maxShift-2))
    
    crawl_PWM = crawl[np.logical_and(crawl["Shift"] >= -minShift + ltrim, crawl['Shift'] <= maxShift - rtrim)].copy()
    crawl_PWM['Seqs'] = [seq[ltrim:-rtrim] for seq in crawl_PWM['Seqs']]

    plotPWM(basename + "_PWM.png", crawl_PWM["Seqs"], crawl_PWM[scorename])
    plotPWM(basename + "_PWM_RC.png", rc(crawl_PWM["Seqs"]), crawl_PWM[scorename])

    # crawl["Seqs"] = crawl["Seqs"].str.replace("-", "")
    # grouping = crawl.groupby("Shift")
    # for group in grouping.groups.keys():
        # temp = grouping.get_group(group)
        # temp.to_csv(basename + "_aligned_shift" + str(group) + ".tsv", index=False, sep="\t")

    summary = crawl['Shift'].value_counts(sort=False).sort_index()
    summary = summary.rename_axis("Shift")
    summary.to_csv(basename + "_TDC_summary.tsv", header=['# Sequences'], sep="\t", index=True)
    
def main():
    TDC(sys.argv[1])
