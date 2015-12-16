from math import sqrt

from PyParameters import Wi, Wf, tau

# cdef double cWi = Wi #3e11
# cdef double cWf = Wf #1e11

N = 1000
g13 = 2.5e9
g24 = 2.5e9
Delta = -2.0e10;
alpha = 1.1e7;
Ng = sqrt(N)* g13;

def Ui(W):
    return -2*(g24**2/Delta)*((Ng**2)*(W**2)/((Ng**2+W**2)**2))

def Jij(Wi, Wj):
    return alpha*Wi*Wj/(sqrt(Ng**2 + Wi**2)* sqrt(Ng**2 + Wj**2))

U0 = Ui(Wi)

mu = 0.5*U0

xi = [1.01015958087519, 0.914144976064563, 1.04162956443615, 1.0679898084607, 0.958180948719382]

U00 = Ui(Wi)

# cdef double ctau = tau #1e-11

def dU0(i):
    return Ui(Wi*xi[i]) - U00

def J0(i, j):
    return Jij(Wi*xi[i], Wi*xi[j])

def Wt(t):
    if t<tau:
        return Wi + (Wf - Wi) * t / tau
    else:
        return Wf + (Wi - Wf) * (t - tau) / tau

def Wtp(t):
    if t<tau:
        return (Wf - Wi) / tau
    else:
        return (Wi - Wf) / tau

# cdef double Wtc(double t):
#     if t<ctau:
#         return cWi + (cWf - cWi) * t / ctau
#     else:
#         return cWf + (cWi - cWf) * (t - ctau) / ctau
#
# cdef double Wtpc(double t):
#     if t<ctau:
#         return (cWf - cWi) / ctau
#     else:
#         return (cWi - cWf) / ctau

def Wit(i, t):
    return Wt(t) * xi[i]

def Witp(i, t):
    return Wtp(t) * xi[i]

def J(i, j, t):
    return Jij(Wit(i, t), Wit(j, t))

def Jp(i, j, t):
    return -((alpha *Wit(i, t)**2 *Wit(j, t) *Witp(i, t))/((Ng**2+Wit(i, t)**2)**(3./2) *sqrt(Ng**2+Wit(j, t)**2)))+(alpha *Wit(j, t) *Witp(i, t))/(sqrt(Ng**2+Wit(i, t)**2)* sqrt(Ng**2+Wit(j, t)**2))-(alpha *Wit(i, t) *Wit(j, t)**2 *Witp(j, t))/(sqrt(Ng**2+Wit(i, t)**2)* (Ng**2+Wit(j, t)**2)**(3./2))+(alpha *Wit(i, t) *Witp(j, t))/(sqrt(Ng**2+Wit(i, t)**2)* sqrt(Ng**2+Wit(j, t)**2))

def U0(t):
    return Ui(Wt(t))

def U0p(t):
    return (8 *g24**2 *Ng**2 *Wt(t)**3 *Wtp(t))/(Delta *(Ng**2+Wt(t)**2)**3)-(4 *g24**2 *Ng**2 *Wt(t) *Wtp(t))/(Delta *(Ng**2+Wt(t)**2)**2)

def dU(i, t):
    return Ui(Wit(i, t)) - U0(t)

# def

# def dU(i, t)