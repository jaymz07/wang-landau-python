import string
import numpy

def fitPolynomial(data2D,degree):
    a, b = [], []
    j=0
    for k in range(0,len(data2D)):
        data=data2D[k]
        ln=[0.0]*(degree+1)
        if(isinstance(data,list)):
            for i in range(0,degree+1):
                ln[i] = float(pow(data[0],i))
        else:
            for i in range(0,degree+1):
                ln[i] = float(pow(k,i))
        y = 0.0
        if(isinstance(data,list)):
            y = float(data[1])
        else:
            y=float(data)
        if(j!=0):
            a[j].append(ln)
            b[j].append(float(y))
        else:
            a.append(ln)
            b.append(float(y))

    
    A, B = numpy.matrix(a), numpy.matrix(b).transpose()
    
    coeff= ((A.transpose() * A).getI()) * A.transpose() * B
    
    return numpy.asarray(coeff)

def fitParabola(data2D):
    s11, s12, s22, sy1, sy2, s1, s2, sy = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    N=len(data2D)
    for i in range(0,N):
        data = data2D[i]
        x, y = 0.0, 0.0
        if(isinstance(data,list)):
            x, y = data[0], data[1]
        else:
            x = i
            y= data
        x1, x2 = x, x**2
        s11 += x1**2
        s12 += x1*x2
        s22 += x2**2
        sy1 += y*x1
        sy2 += y*x2
        s1 += x1
        s2 += x2
        sy += y
    s11 -= s1**2/N
    s12 -= s1*s2/N
    s22 -= s2**2/N
    sy1 -= sy*s1/N
    sy2 -= sy*s2/N
    
    out2 = (sy1*s22-sy2*s12)/(s22*s11-s12**2)
    out3 = (sy2*s11 - sy1* s12)/(s22*s11 - s12**2)
    out1 = sy/N - out2*s1/N - out3*s2/N
    
    return [out1,out2, out3]