N = 1000
a = 1
G = 6.67 * 10E-11
M = 2*10E30

h = a/(N-1)

t = [k*h for k in range(1001)]


# Vx = x'

# Vx' = - GM(Ax(t))  

# Vy = y'

# Vy' =  - GM(Ay(t)) 

zx = [0]
zy = [11.5]

x = [0.5]

y = [0] 



def fx(tk) : 
    return - G*M*x[tk]/((x[tk]**2 + y[tk]**2)**3/2)

def fy(tk) : 
    return - G*M*y[tk]/((x[tk]**2 + y[tk]**2)**3/2)


for tk in range(N) :
    x.append(x[-1]+h*zx[-1]) 

    zx.append(zx[-1] + h*fx(tk))

    y.append(y[-1]+h*zy[-1]) 

    zy.append(zy[-1] + h*fy(tk))


from matplotlib import pyplot as plt

plt.plot(x,y)



