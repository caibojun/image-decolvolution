from astropy.io import fits
import numpy as np

########################################################################

def ProjFletNon(b,c,dia):
	Lambda=0
	dlambda=1
	tol_r=10**(-11)*b
	tol_lam=10**(-11)
	biter=0
	siter=0
	maxprojections=1000
	
	x=np.array([map(lambda m,n:max(m,n),m,n) for (m,n) in zip(np.ones([c.shape[0],1]),(c+Lambda)/dia)])
	r=np.sum(x)-b
	if abs(r)<tol_r:
		return c/dia
	if r<0:
		lambdal=Lambda
		rl=r
		Lambda+=dlambda
		x=np.array([map(lambda m,n:max(m,n),m,n) for (m,n) in zip(np.zeros([c.shape[0],1]),(c+Lambda)/dia)])
		r=np.sum(x)-b
		while r<0:
			biter+=1
			lambdal=Lambda
			s=max(rl/r-1,0.1)
			dlambda+=dlambda/s
			Lambda+=dlambda
			rl=r
			x=np.array([map(lambda m,n:max(m,n),m,n) for (m,n) in zip(np.zeros([c.shape[0],1]),(c+Lambda)/dia)])
			r=np.sum(x)-b
		lambdau=Lambda
		ru=r
	else:
		lambdau=Lambda
		ru=r
		Lambda-=dlambda
		x=np.array([map(lambda m,n:max(m,n),m,n) for (m,n) in zip(np.zeros([c.shape[0],1]),(c+Lambda)/dia)])
		r=np.sum(x)-b
		while r>0:
			biter+=1
			lambdau=Lambda
			s=max(ru/r-1,0.1)
			dlambda+=dlambda/s
			Lambda-=dlambda
			ru=r
			x=np.array([map(lambda m,n:max(m,n),m,n) for (m,n) in zip(np.zeros([c.shape[0],1]),(c+Lambda)/dia)])
			r=sum(x)-b
		lambdal=Lambda
		rl=r
	if abs(ru)<tol_r:
		x=np.array([map(lambda m,n:max(m,n),m,n) for (m,n) in zip(np.zeros([c.shape[0],1]),(c+lambdau)/dia)])
	if abs(rl)<tol_r:
		x=np.array([map(lambda m,n:max(m,n),m,n) for (m,n) in zip(np.zeros([c.shape[0],1]),(c+lambdal)/dia)])
	s=1-rl/ru
	dlambda/=s
	Lambda=lambdau-dlambda
	x=np.array([map(lambda m,n:max(m,n),m,n) for (m,n) in zip(np.zeros([c.shape[0],1]),(c+Lambda)/dia)])
	r=np.sum(x)-b
	maxit_s=maxprojections-biter
	
	while abs(r)>tol_r and(dlambda>tol_lam*(1+abs(Lambda)))and siter<maxit_s:
		siter+=1
		if r>0:
			if s<2:
				lambdau=Lambda
				ru=r
				s=1-rl/ru
				dlambda=(lambdau-lambdal)/s
				Lambda=lambdau-dlambda
			else:
				s=max(ru/r-1,0.1)
				dlambda=(lambdau-lambdal)/s
				lambda_new=max(Lambda-dlambda,0.75*lambdal+0.25*Lambda)
				lambdau=Lambda
				ru=r
				Lambda=lambda_new
				s=(lambdau-lambdal)/(lambdau-Lambda)
		else:
			if s>2:
				lambdal=Lambda
				rl=r
				s=1-rl/ru
				dlambda=(lambdau-lambdal)/s
				Lambda=lambdau-dlambda
			else:
				s=max(rl/r-1,0.1)
				dlambda=(Lambda-lambdal)/s
				lambda_new=min(Lambda-dlambda,0.75*lambdal+0.25*Lambda)
				lambdal=Lambda
				rl=r
				Lambda=lambda_new
				s=(lambdau-lambdal)/(lambdau-Lambda)
		x=np.array([map(lambda m,n:max(m,n),m,n) for (m,n) in zip(np.zeros([c.shape[0],1]),(c+Lambda)/dia)])
		r=np.sum(x)-b
	return x

########################################################################

def fourierfun(X,TF):
	Xpd=np.fft.fft2(X)
	Xoverpad=np.fft.ifft2(Xpd*TF)
	X=Xoverpad.real
	return X

########################################################################	

def sgpdecon(paralist):
	"""Input parameters
	"""
	originalimage=paralist[0]
	originalpsf=paralist[1]
	iteration=paralist[2]
	"""Default
	"""
	alpha_min=10**(-5)
	alpha_max=10**5
	zeta=10**(-4)
	initalpha=1.3
	M=1
	Malpha=3
	tau=0.5
	"""Change the shape of psf
	"""
	Inpsf=originalpsf
	TF = np.fft.fft2(np.fft.fftshift(Inpsf))
	TF = TF.reshape(len(TF)**2,-1)
	CTF = np.conj(TF)
	"""Change the shape of fits img
	"""
	x=originalimage.reshape(len(originalimage)**2,-1)
	gn=x.copy()
	flux=np.sum(x)
	Iter=1
	alpha=initalpha
	Valpha=alpha_max*np.ones([Malpha,1])
	Fold=np.ones([M,1])
	ONE=np.ones([gn.shape[0],1])
	x=ProjFletNon(flux,x,np.ones([x.shape[0],1]))
	x_tf=fourierfun(x,TF)
	den=x_tf
	temp=gn/den
	g=ONE-fourierfun(temp,CTF)
	fv=np.sum(gn*np.log(temp))+np.sum(x_tf)-flux
	y=fourierfun(gn,CTF)
	x_low_bound=min(y)
	x_up_bound=max(y)
	if x_up_bound/x_low_bound<50:
		x_low_bound/=10
		x_up_bound*=10
	xx=x
	xx=np.array([map(lambda m,n:max(m,n),m,n )for (m,n) in zip(xx,x_low_bound*np.ones([xx.shape[0],1]))])
	xx=np.array([map(lambda m,n:min(m,n),m,n )for (m,n) in zip(xx, x_up_bound*np.ones([xx.shape[0],1]))])
	D=1/xx
	loop=1
	while loop==1:
		Valpha[0:Malpha-3]=Valpha[1:Malpha-2]
		if M>1:
			Fold[0:M-3]=Fold[1:M-2]
		Fold[M-1]=fv
		y=x-alpha*xx*g
		y=ProjFletNon(flux,y*D,D)
		d=y-x
		gd=np.sum(d*g)
		lam=1
		fcontinue=1
		d_tf=fourierfun(d,TF)
		fr=np.max(Fold)
		while fcontinue==1:
			xplux=x+lam*d
			x_tf_try=x_tf+lam*d_tf
			den=x_tf_try+0
			temp=gn/den
			# Estimate whether the probability distribution in KL divergence is positive
			if np.sum(temp)==np.sum(abs(temp)):
				fv=np.sum(gn*np.log(temp))+np.sum(x_tf)-flux
				if lam<10**(-12) or fv<=(fr+zeta*lam*gd):
					x=xplux
					print x
					sk=lam*d
					gtemp=ONE-fourierfun(temp,CTF)
					yk=gtemp-g
					g=gtemp
					fcontinue=0
				else:
					lam*=zeta
			else:
				lam*=zeta
		xx=x
		xx=np.array([map(lambda m,n:max(m,n),m,n )for (m,n) in zip(xx,x_low_bound*np.ones([xx.shape[0],1]))])
		xx=np.array([map(lambda m,n:min(m,n),m,n )for (m,n) in zip(xx, x_up_bound*np.ones([xx.shape[0],1]))])
		D=np.ones([xx.shape[0],1])/xx
		sk2=sk*D
		yk2=yk*xx
		bk=sum(sk2*yk)
		ck=sum(yk2*sk)
		if bk<=0:
			alpha1=min(10*alpha,alpha_max)
		else:
			alpha1BB=sum(sk2*sk2)/bk
			alpha1=min(alpha_max,max(alpha_min,alpha1BB))
		if ck<=0:
			alpha2=min(10*alpha,alpha_max)
		else:
			alpha2BB=ck/sum(yk2*yk2)
			alpha2=min(alpha_max,max(alpha_min,alpha2BB))
		if Iter<=20:
			alpha=np.min(Valpha)
		else:
			if alpha2/alpha1<tau:
				alpha=np.min(Valpha)
				tau*=0.9
			else:
				alpha=alpha1
				tau*=1.1
		Iter+=1
		if Iter>iteration:
			loop=0
		print "x=",x
	restoredimage=x.reshape(originalimage.shape[0],originalimage.shape[1])
	return restoredimage
