"""for Python"""
import sys,os,glob ,math,scipy,itertools,time
import numpy as np
from astropy.io import fits
#import matplotlib.pyplot as plt

####################bwlabel.py#####################
"""this part is disign for deblending region from pics have been seperated background already"""

def size(IN):
        M = len(IN[0])
        N = len(IN)
        return (M, N)

def NumberOfRuns(IN):
        M, N = size(IN)
        result = 0
        if M != 0 and N != 0:
                for col in IN:
                        if col[0] != 0:
                                result += 1
                        for idx in range(1, M):
                                if col[idx] != 0 and col[idx-1] == 0:
                                        result += 1
        return result

def FillRunVectors(IN):
        M, N = size(IN)
        c = []
        sr = []
        er = []
        for cidx ,col in enumerate(IN):
                k = 0
                while k < M:
                        try:
                                k += col[k:].index(1)
                                c.append(cidx+1)
                                sr.append(k+1) #! for matlab
                                try:
                                        k += col[k:].index(0)
                                except ValueError:
                                        k = M
                                er.append(k)
                        except ValueError:
                                break
        return sr, er, c
           
def FirstPass(numRuns, mode, sr, er, c):
        currentColumn = 0
        nextLabel = 1
        firstRunOnPreviousColumn = -1
        lastRunOnPreviousColumn = -1
        firstRunOnThisColumn = -1
        equivList = []
        labels = [0] * numRuns
        if mode == 8:
                offset = 1
        else:
                offset = 0
        for k in range(numRuns):
                if c[k] == currentColumn + 1:
                        firstRunOnPreviousColumn = firstRunOnThisColumn
                        firstRunOnThisColumn = k
                        lastRunOnPreviousColumn = k-1
                        currentColumn = c[k]
                elif c[k] > currentColumn+1:
                        firstRunOnPreviousColumn = -1
                        lastRunOnPreviousColumn = -1
                        firstRunOnThisColumn = k
                        currentColumn = c[k]
                else:
                        pass
                if firstRunOnPreviousColumn >= 0:
                        p = firstRunOnPreviousColumn
                        while p<=lastRunOnPreviousColumn and sr[p]<=er[k]+offset:
                                if er[k]>=sr[p]-offset and sr[k]<=er[p]+offset:
                                        if labels[k] == 0:
                                                labels[k] = labels[p]
                                        else:
                                                if labels[k] != labels[p]:
                                                        equivList.insert(0, (labels[k],labels[p]))
                                                else:
                                                        pass
                                p += 1
                if labels[k] == 0:
                        labels[k] = nextLabel
                        nextLabel += 1
        rowEquivalences = []
        colEquivalences = []
        if len(equivList) > 0:       
                for item0, item1 in equivList:
                        rowEquivalences.append(item0)
                        colEquivalences.append(item1)
        return labels, rowEquivalences, colEquivalences

def bwlabel1(BW, mode=8):
        BW = list(BW)   
        for k in range(len(BW)):
                BW[k] = list(BW[k])
        numRuns = NumberOfRuns(BW)
        sr, er, c = FillRunVectors(BW)
        labels, rowEquivalences, colEquivalences = FirstPass(numRuns, mode, sr, er, c)
        return sr, er, c, labels, rowEquivalences, colEquivalences

def bwlabel1p5(labels, rowEq, colEq):
        lblist = []
        for k, r in enumerate(rowEq):
                if r==-1:
                        continue
                queue = list([rowEq[k], colEq[k]])
                cur_lblist = list()
                while queue:
                        head = queue.pop(0)
                        cur_lblist.append(head)
                        for n in range(k+1, len(rowEq)):
                                if rowEq[n] == head:
                                        queue.append(colEq[n])
                                        rowEq[n] = -1
                                        colEq[n] = -1
                                elif colEq[n] == head:
                                        queue.append(rowEq[n])
                                        rowEq[n] = -1
                                        colEq[n] = -1
                lblist.append(cur_lblist)
# for test
#    print lblist     # print wrong label
        labels = scipy.array(labels)
        for oldlabels in lblist:
                for ol in oldlabels:
                        labels[labels==ol] = oldlabels[0]
        sort_labels = zip(labels, range(len(labels)))
        sort_labels.sort()
        sort_idx = [k[1] for k in sort_labels]
        sort_labels = scipy.array([k[0] for k in sort_labels])
        if sort_labels[0] != 1:
                sort_labels -= (sort_labels[0]-1)
        for k in range(1, len(sort_labels)):
                cur_label = sort_labels[k]
                pre_label = sort_labels[k-1]
                if cur_label > pre_label+1:
                        sort_labels[sort_labels==cur_label] = pre_label + 1
        for k, l in itertools.izip(sort_idx, sort_labels):
                        labels[k] = l
#    print labels
        return labels

def bwlabel2(sr, er, sc, labels, M, N):
        L = scipy.zeros((M, N))
        for s, e, c, l in itertools.izip(sr, er, sc, labels):
                L[c-1, (s-1):e] = l
        return L

def bwlabel(BW, mode=8):
        sr, er, sc, labels, req, ceq = bwlabel1(BW, mode)
        labels = bwlabel1p5(labels, req, ceq)
        M, N = size(BW)
        return bwlabel2(sr, er, sc, labels, N, M)
############################################################

def Open_files():
        """Indicate that "Input the FITS Path:",then put the gray level imfornation in to matrix.

        Parameter
        ---------
        Path: a string variable is used to input file path.
        Gray_level: ndarray , a list matrix variable store the gray value.
    
        P.S.
        when you cannot input a correct files path or you want to quit the program , input "exit"
        then the program would shutdown."""
    
        

        Path = raw_input("Input file path:")
        Path1 = Path + '/*.fit*'
        print Path1
        if __name__ == '__main__':
                fitfile = glob.glob(Path1)
        while len(fitfile)==0:
                print "Error:File doesn't exist,please check the path.\n"
                Path = raw_input("Input correct path:")
                if Path == "exit":
                        return
                elif Path == '':
                        Path=sys.path[0]
                Path1 = Path + '/*.fit*'
                if __name__ == '__main__':
                        fitfile = glob.glob(Path1)
        Save = Path[:len(Path)-len(Path.split("/")[-1])] + 'star'   
        if os.path.isdir(Save)==True:
                os.popen('rm -rf '+Save)
                os.mkdir(Save)
        elif Save=='':
                Save=sys.path[0]
        else:
                os.mkdir(Save)
        mesh=raw_input("Input the mesh size :")
        fitfile.sort()
        for couts,fit in enumerate(fitfile):
                print couts+1,fit
#               print "----------------------------------------\n"
                hdu = fits.open(fit)
                gray = hdu[0].data
                h=hdu[0].header['NAXIS1']
                w=hdu[0].header['NAXIS2']
                if h!=w:
					h=min(h,w)
					w=h
                gray = np.transpose(gray)
                s1=time.clock()
                starlocat=[]
                list4=Mesh(Save,fit,mesh,gray,h,w,starlocat)
                stararray=np.array(list4,dtype=int)
                np.savetxt(Save+'/'+fit.split("/")[-1]+'.txt',stararray,fmt="%d",delimiter="\t")
                e1=time.clock()
                print "Single picture spends: %f seconds"%(e1-s1)

def Mesh(Save,fit,window,gray,h,w,list1):
        """create a background map ,then subtracted from imagine.

        Parameter
        ---------
        mesh: set a value as a window to cut the whole pic, the value we selected must can divide the "h" and "w". 
        mesh_w: X aixs boundary of the mesh.
        mesh_h: Y aixs boundary of the mesh.
        threshold: as a extractor.for 7802 data we choose 2.0
        """
        mesh_w=0
        mesh_h=0
        i=0
        j=0
        copy=gray.copy()
        mesh=int(window)
        mesh_w=i+mesh
        mesh_h=j+mesh
        while mesh_h<=h:
                while mesh_w<=w:
                        narry = gray[i:mesh_w,j:mesh_h]
                        mean=np.mean(narry)
                        median=np.median(narry)
                        std=np.std(narry)
                        threshold=median+3.0*std
                        narry[narry<threshold]=0
                        narry[narry>=threshold]=1
                        gray[i:mesh_w,j:mesh_h]=narry
                        i+=mesh
                        mesh_w+=mesh
                i=0
                mesh_w=i+mesh
                j+=mesh
                mesh_h+=mesh
        img = bwlabel(gray)
        list3=pickout(Save,fit,img,copy,list1)
        return list3

def pickout(Save,fit,mat,orimat,list1):
        #filepath=fit.split("/")[-1]
        #f=open(Save+'/'+filepath+'.txt','w+')
        num=int(np.amax(mat))
        for i in range(1,num+1):
                count=np.sum(mat==i)
                obj=mat.copy()
                obj[mat==i]=1
                obj[mat!=i]=0
                obj=orimat*obj
                count=np.sum(obj!=0)
                if count>=2:#delete regions too small to extract PSFs
                        new_thresh=0.5*np.max(obj)
                        obj[obj<=new_thresh]=0
                        obj[obj>new_thresh]=1
                        list2=pickoutv1(Save,fit,obj,orimat,list1)
        return list2

def pickoutv1(Save,fit,mat,orimat,list1):
        star = bwlabel(mat)
        labels=1
        xval=0
        yval=0
        number=int(np.amax(star))
        for i in range(1,number+1):
                clusters=star.copy()
                clusters[clusters==i]=1
                clusters[clusters!=i]=0
                clusters*=orimat
                counter=np.sum(clusters!=0)
#                if counter>=2:#delete regions too small to extract PSFs
                element=np.where(clusters==np.max(clusters))
                xval=element[0][0]+1
                yval=element[1][0]+1
                location=[xval,yval]
                list1.append(location)
        return list1
Open_files()
