import numpy as np
from matplotlib import pyplot as plt
import copy

# defining properties of fluid and dimensions of the domain
rho = 1 ; mu = 2*10**-5; L = 0.2 ; l = 0.01 ; h1 = 0.01 ; exprat = 2 ; h2 = exprat*h1 ; dx =0.005 ; dy =0.005 ; ax = 0 ; ay = -9.81; G = 0
alpha_p = 0.5 ; alpha_vel = 0.5;alpha_prime =1;alpha = 0.5
#uniform structured grid creation
nx = int(L/dx) ;ny = int(h2/dy) ; ny1 = int((h2-h1)/dy) ; nx1 = int((l/dx))
x = np.linspace(0,L+2*dx,nx+2)
y = np.linspace(0,h2+2*dy,ny+2)
X,Y = np.meshgrid(x,y)
# inlet velocity (parabolic)
y1 = np.linspace(dy/2, h1, ny-ny1)
Ui = 2*(y1/h1)*(1-(y1/h1))
#guess values
p = np.zeros((ny + 2, nx + 2)) ; U = np.zeros((ny + 2, nx + 2)) ; V = np.zeros((ny + 2, nx + 2)) ; Uold = np.zeros((ny + 2, nx + 2)) ; Vold = np.zeros((ny + 2, nx + 2)) ; pold = np.zeros((ny + 2, nx + 2))
Uf = np.ones((ny+2,nx+3))/10 ; Vf = np.zeros((ny+3,nx+2));
#initializing link coefficents
ae =np.zeros((ny + 2, nx + 2));ap =np.zeros((ny + 2, nx + 2));aw =np.zeros((ny + 2, nx + 2));an =np.zeros((ny + 2, nx + 2));aS =np.zeros((ny + 2, nx + 2));Sx = np.zeros((ny + 2, nx + 2));Sy = np.zeros((ny + 2, nx + 2))
Fx =rho*Uf*dy ; Fy = rho*Vf*dx ;  Dx = (mu*dy)/dx ; Dy = (mu*dx)/dy;
p_prime = np.zeros((ny + 2, nx + 2));p_primeold = np.zeros((ny + 2, nx + 2))
def link_coefficients():

    # interior cells block-1
    for i in range(ny1+2,ny):
        for j in range(2,nx1+1):
            ae[i][j] = Dx + max(-Fx[i][j+1],0)
            aw[i][j] = Dx + max(Fx[i][j],0)
            an[i][j] = Dy + max(-Fy[i+1][j],0)
            aS[i][j] = Dy + max(Fy[i][j],0)
            ap[i][j] = 2*Dx+2*Dy +max(Fx[i][j+1],0)+max(-Fx[i][j],0)+max(-Fy[i][j],0)+max(Fy[i+1][j],0)
            Sx[i][j] = (p[i][j - 1] - p[i][j + 1]) * dy * 0.5
            Sy[i][j] = (p[i - 1][j] - p[i + 1][j]) * dx * 0.5 + G
    # interior block-2
    for i in range(2,ny):
        for j in range(nx1+2,nx):
            ae[i][j] = Dx + max(-Fx[i][j + 1], 0)
            aw[i][j] = Dx + max(Fx[i][j], 0)
            an[i][j] = Dy + max(-Fy[i + 1][j], 0)
            aS[i][j] = Dy + max(Fy[i][j], 0)
            ap[i][j] = 2 * Dx + 2 * Dy + max(Fx[i][j + 1], 0) + max(-Fx[i][j], 0) + max(-Fy[i][j], 0) + max(Fy[i + 1][j], 0)
            Sx[i][j] = (p[i][j - 1] - p[i][j + 1]) * dy * 0.5
            Sy[i][j] = (p[i - 1][j] - p[i + 1][j]) * dx * 0.5 + G
    #  inlet boundary except corners
    for i in range(ny1+2,ny):
        j = 1
        ae[i][j] = Dx + max(0.0, -Fx[i][j+1])
        aw[i][j] = 0
        an[i][j] = Dy + max(0.0, -Fy[i+1][j])
        aS[i][j] = Dy + max(0.0, Fy[i][j])
        ap[i][j] = 3*Dx+2*Dy + max(0.0, Fx[i][j+1]) + max(0.0, -Fx[i][j]) + max(0.0, Fy[i+1][j]) + max(0.0, -Fy[i][j])
        Sx[i][j] = 0.5 * (p[i][j] - p[i][j+1]) * dy + rho * Uf[i][j] * Uf[i][j] * dy + 2 * Dx * Uf[i][j]
        Sy[i][j] = (p[i - 1][j] - p[i + 1][j]) * dx * 0.5 + G
    # top wall
    i=ny
    for j in range(2,nx):
        ae[i][j] = Dx + max(-Fx[i][j + 1], 0)
        aw[i][j] = Dx + max(Fx[i][j], 0)
        an[i][j] = 0
        aS[i][j] = Dy + max(0.0, Fy[i][j])
        ap[i][j] = Dx + Dx + Dy + max(0.0, Fx[i][j+1]) + max(0.0, -Fx[i][j]) + max(0.0, Fy[i+1][j]) + max(0.0, -Fy[i][j]) + 2 * Dy
        Sx[i][j] = (p[i][j - 1] - p[i][j + 1]) * dy * 0.5
        Sy[i][j] = (p[i - 1][j] - p[i][j]) * dx * 0.5 + G
    # bottom wall in first block
    i = ny1+1
    for j in range(2,nx1+1):
        ae[i][j] = Dx + max(-Fx[i][j + 1], 0)
        aw[i][j] = Dx + max(Fx[i][j], 0)
        an[i][j] = Dy + max(-Fy[i + 1][j], 0)
        aS[i][j] = 0
        ap[i][j] = Dx + Dx + Dy + max(0.0, Fx[i][j+1]) + max(0.0, -Fx[i][j]) + max(0.0, Fy[i+1][j]) + max(0.0, -Fy[i][j]) + 2 * Dy
        Sx[i][j] = (p[i][j - 1] - p[i][j + 1]) * dy * 0.5
        Sy[i][j] = (p[i][j] - p[i + 1][j]) * dx * 0.5 + G
    # wall on the step
    j = nx1+1
    for i in range(2,ny1+1):
        ae[i][j] = Dx + max(-Fx[i][j + 1], 0)
        aw[i][j] = 0
        an[i][j] = Dy + max(-Fy[i + 1][j], 0)
        aS[i][j] = Dy + max(0.0, Fy[i-1][j])
        ap[i][j] = Dx + Dy + Dy + max(0.0, Fx[i][j+1]) + max(0.0, -Fx[i][j-1]) + max(0.0, Fy[i+1][j]) + max(0.0, -Fy[i-1][j]) + 2 * Dx
        Sx[i][j] = (p[i][j] - p[i][j + 1]) * dy * 0.5
        Sy[i][j] = (p[i - 1][j] - p[i + 1][j]) * dx * 0.5 + G
    # outlet_face
    j= nx+1
    for i in range(2,ny):
        ae[i][j] = 0
        aw[i][j] = Dx + max(Fx[i][j], 0)
        an[i][j] = Dy + max(-Fy[i + 1][j], 0)
        aS[i][j] = Dy + max(0.0, Fy[i - 1][j])
        ap[i][j] = Dx + Dy + Dy + max(0.0, Fx[i][j + 1]) + max(0.0, -Fx[i][j - 1]) + max(0.0, Fy[i + 1][j]) + max(0.0, -Fy[i - 1][j]) + 2 * Dx
        Sx[i][j] = (p[i][j-1] + p[i][j ]) * dy * 0.5
        Sy[i][j] = (p[i - 1][j] - p[i + 1][j]) * dx * 0.5 + G

    # bottom wall in block-2
    i=1
    for j in range(nx1+2,nx):
        ae[i][j] = Dx + max(-Fx[i][j + 1], 0)
        aw[i][j] = Dx + max(Fx[i][j], 0)
        an[i][j] = Dy + max(-Fy[i + 1][j], 0)
        aS[i][j] = 0
        ap[i][j] = Dx + Dx + Dy + max(0.0, Fx[i][j+1]) + max(0.0, -Fx[i][j]) + max(0.0, Fy[i+1][j]) + max(0.0, -Fy[i][j]) + 2 * Dy
        Sx[i][j] = 0.5 * (p[i][j-1] - p[i][j+1]) * dy
        Sy[i][j] = 0.5 * (p[i][j] - p[i+1][j]) * dx +G
    j =nx
    for i in range(2,ny):
        ae[i][j] = 0
        aw[i][j] = Dx + max(Fx[i][j], 0)
        an[i][j] = Dy + max(-Fy[i + 1][j], 0)
        aS[i][j] = Dy + max(0.0, Fy[i][j])
        ap[i][j] = Dx + Dy + Dy + max(0.0, Fx[i][j+1]) + max(0.0, -Fx[i][j]) + max(0.0, Fy[i+1][j]) + max(0.0, -Fy[i][j])
        Sx[i][j] = 0.5 * (p[i][j] + p[i][j-1]) * dy
        Sy[i][j] = 0.5 * (p[i-1][j] - p[i+1][j]) * dx +G
    #connecting 2 blocks cells
    j = nx1+1
    for i in range(ny1+1,ny):
        ae[i][j] = Dx + max(-Fx[i][j + 1], 0)
        aw[i][j] = Dx + max(Fx[i][j ], 0)
        an[i][j] = Dy + max(-Fy[i + 1][j], 0)
        aS[i][j] = Dy + max(Fy[i ][j], 0)
        ap[i][j] = Dx + Dx + Dy + Dy + max(0.0, Fx[i][j+1]) + max(0.0, -Fx[i][j]) + max(0.0, Fy[i+1][j]) + max(0.0, -Fy[i][j])
        Sx[i][j] = 0.5 * (p[i][j-1] - p[i][j+1]) * dy
        Sy[i][j] = 0.5 * (p[i-1][j] - p[i+1][j]) * dx
    # top corner
    j = 1
    i = ny
    ae[i][j] = Dx + max(-Fx[i][j + 1], 0)
    aw[i][j] = 0
    an[i][j] = 0
    aS[i][j] = Dy + max(Fy[i][j], 0)
    ap[i][j] = Dx + Dy + max(0.0, Fx[i][j+1]) + max(0.0, -Fx[i][j]) + max(0.0, Fy[i+1][j]) + max(0.0, -Fy[i][j]) + 2 * Dy + 2 * Dx
    Sx[i][j] = 0.5 * (p[i][j] - p[i][j+1]) * dy + rho * Uf[i][j] * Uf[i][j]*dy + 2 * Dx * Uf[i][j]
    Sy[i][j] = 0.5 * (p[i-1][j] - p[i][j]) * dx +G
    # bottom corner block-1
    j = 1
    i = ny1+1
    ae[i][j] = Dx + max(-Fx[i][j + 1], 0)
    aw[i][j] = 0
    an[i][j] = Dy + max(-Fy[i + 1][j], 0)
    aS[i][j] = 0
    ap[i][j] = Dx + Dy + max(0.0, Fx[i][j+1]) + max(0.0, -Fx[i][j]) + max(0.0, Fy[i+1][j]) + max(0.0, -Fy[i][j]) + 2 * Dy + 2 * Dx

    Sx[i][j] = 0.5 * (p[i][j] - p[i][j+1]) * dy + rho * Uf[i][j] * Uf[i][j] *dy + 2 * Dx * Uf[i][j]
    Sy[i][j] = 0.5 * (p[i][j] - p[i+1][j]) * dx + G

    # top right corner
    j = nx
    i = ny

    ae[i][j] = 0
    aw[i][j] = Dx + max(Fx[i][j], 0)
    an[i][j] = 0
    aS[i][j] = Dy + max(Fy[i][j], 0)
    ap[i][j] = Dx + Dy +  max(0.0, Fx[i][j+1]) + max(0.0, -Fx[i][j]) + max(0.0, Fy[i+1][j]) + max(0.0, -Fy[i][j]) + 2 * Dy
    Sx[i][j] = 0.5 * (p[i][j] + p[i][j-1]) * dy
    Sy[i][j] = 0.5 * (p[i-1][j] - p[i][j]) * dx +G

    # bottom right corner
    j = nx;
    i = 1;


    ae[i][j] = 0;
    aw[i][j] = Dx + max(Fx[i][j], 0)
    an[i][j] = Dy + max(-Fy[i + 1][j], 0)
    aS[i][j] = 0;
    ap[i][j] = Dx + Dy + max(0.0, Fx[i][j+1]) + max(0.0, -Fx[i][j]) + max(0.0, Fy[i+1][j]) + max(0.0, -Fy[i][j]) + 2 * Dy;
    Sx[i][j] = 0.5 * (p[i][j] + p[i][j-1]) * dy;
    Sy[i][j] = 0.5 * (p[i][j] - p[i+1][j]) * dx +G

    # corner cell near step
    j = nx1 + 1
    i = 1

    ae[i][j] = Dx + max(-Fx[i][j + 1], 0)
    aw[i][j] = 0
    an[i][j] = Dy + max(-Fy[i + 1][j], 0)
    aS[i][j] = 0
    ap[i][j] = Dx + Dy + max(0.0, Fx[i][j+1]) + max(0.0, -Fx[i][j]) + max(0.0, Fy[i+1][j]) + max(0.0, -Fy[i][j]) + 2 * Dy + 2 * Dx
    Sx[i][j] = 0.5 * (p[i][j] - p[i][j+1]) * dy
    Sy[i][j] = 0.5 * (p[i][j] - p[i+1][j]) * dx +G
def solve_xandy_momentum():
    # solving x and y momentum equations
    Uold = copy.deepcopy(U)
    Vold = copy.deepcopy(V)
    for k in range(1,500):
        for i in range(1,ny1+1):
            for j in range(nx1+1,nx+1):
                U[i][j] = alpha*((ae[i][j]*U[i][j+1] + aw[i][j]*U[i][j-1]+aS[i][j]*U[i-1][j]+an[i][j]*U[i+1][j]+Sx[i][j])/ap[i][j])+(1-alpha)*U[i][j]
                V[i][j] = alpha*((ae[i][j] * V[i][j + 1] + aw[i][j] * V[i][j - 1] + aS[i][j] * V[i - 1][j] + an[i][j] * V[i + 1][j] +Sy[i][j])) / ap[i][j] +(1-alpha)*U[i][j]
        for i in range(ny1+1,ny+1):
            for j in range(1,nx+1):
                U[i][j] = alpha*(ae[i][j] * U[i][j + 1] + aw[i][j] * U[i][j - 1] + aS[i][j] * U[i - 1][j] + an[i][j] * U[i + 1][j] + Sx[i][j]) / ap[i][j]+(1-alpha)*U[i][j]
                V[i][j] = alpha*(ae[i][j] * V[i][j + 1] + aw[i][j] * V[i][j - 1] + aS[i][j] * V[i - 1][j] + an[i][j] * V[i + 1][j] + Sy[i][j]) / ap[i][j]+(1-alpha)*V[i][j]
    Res_u = np.linalg.norm(abs(U-Uold))
    Res_v = np.linalg.norm(abs(V-Vold))
    return [Res_u,Res_v]

def face_velocities():
    # velocities on step face
    Uf[1:ny1, nx1 + 1] = 0
    # face velocities
    for j in range(nx1+3,nx):
        for i in range(1,ny1+1):
            Uf[i][j] = 0.5*(U[i][j]+U[i][j-1]) + 0.25*dy*((p[i][j]-p[i][j-2])/ap[i][j-1] + (p[i][j+1]-p[i][j-1])/ap[i][j]) - 0.5*dy*(1/ap[i][j]+1/ap[i][j-1])*(p[i][j]-p[i-1][j])
    for j in range(3,nx):
        for i in range(ny1+1,ny+1):
            Uf[i][j] = 0.5*(U[i][j]+U[i][j-1]) + 0.25*dy*((p[i][j]-p[i][j-2])/ap[i][j-1] + (p[i][j+1]-p[i][j-1])/ap[i][j]) - 0.5*dy*(1/ap[i][j]+1/ap[i][j-1])*(p[i][j]-p[i][j-1])
    #next to inlet
    j = 2;
    for i in range(ny1+1, ny+1):
        Uf[i][j] = 0.5*(U[i][j]+U[i][j-1]) + 0.25*dy*((p[i][j]-p[i][j-1])/ap[i][j-1] + (p[i][j+1]-p[i][j-1])/ap[i][j]) - 0.5*dy*(1/ap[i][j]+1/ap[i][j-1])*(p[i][j]-p[i-1][j])
    # next to step
    j = nx1 + 2;
    for i in range(1, ny1+1):
        Uf[i][j] = 0.5*(U[i][j]+U[i][j-1]) + 0.25*dy*((p[i][j]-p[i][j-1])/ap[i][j-1] + (p[i][j+1]-p[i][j-1])/ap[i][j]) - 0.5*dy*(1/ap[i][j]+1/ap[i][j-1])*(p[i][j]-p[i][j-1])
    #adjacent to outlet
    j = nx
    for i in range(1,ny+1):
        Uf[i][j] = 0.5*(U[i][j]+U[i][j-1]) + 0.25*dy*((p[i][j]-p[i][j-2])/ap[i][j-1] + (p[i][j]-p[i][j-1])/ap[i][j]) - 0.5*dy*(1/ap[i][j]+1/ap[i][j-1])*(p[i][j]-p[i][j-1])

    #v face velocities
    #interior cell faces
    for j in range(1, nx1+1):
        for i in range(ny1+3,ny):
            Vf[i][j] = 0.5*(V[i][j]+V[i-1][j]) + 0.25*dx*((p[i][j]-p[i-2][j])/ap[i-1][j] + (p[i+1][j]-p[i-1][j])/ap[i][j]) - 0.5*dx*(1/ap[i][j]+1/ap[i-1][j])*(p[i][j]-p[i-1][j])
    for j in range (nx1+1, nx+1):
        for i in range (3,ny):
            Vf[i][j] = 0.5 * (V[i][j]+V[i-1][j]) + 0.25 * dx * ((p[i][j]-p[i-1][j]) / ap[i-1][j] + (p[i+1][j]-p[i-1][j]) / ap[i][j]) - 0.5 * dx * (1 / ap[i][j]+1 / ap[i-1][j]) * (p[i][j]-p[i-1][j])
    #faces below top boundary
    i = ny
    for j in range(1,nx+1):
        Vf[i][j] = 0.5*(V[i][j]+V[i-1][j]) + 0.25*dx*((p[i][j]-p[i-2][j])/ap[i-1][j] + (p[i][j]-p[i-1][j])/ap[i][j]) - 0.5*dx*(1/ap[i][j]+1/ap[i-1][j])*(p[i][j]-p[i-1][j])

    #faces above bottom wall block 1
    i= ny1+ 2
    for j in range(1,nx1+1):
        Vf[i][j] = 0.5*(V[i][j]+V[i-1][j]) + 0.25*dx*((p[i][j]-p[i-1][j])/ap[i-1][j] + (p[i+1][j]-p[i-1][j])/ap[i][j]) - 0.5*dx*(1/ap[i][j]+1/ap[i-1][j])*(p[i][j]-p[i-1][j])

    #faces above Bottom wall block 2
    i = 2
    for j in range(nx1+1,nx+1):
        Vf[i][j] = 0.5*(V[i][j]+V[i-1][j]) + 0.25*dx*((p[i][j]-p[i-1][j])/ap[i-1][j] + (p[i+1][j]-p[i-1][j])/ap[i][j]) - 0.5*dx*(1/ap[i][j]+1/ap[i-1][j])*(p[i][j]-p[i-1][j])
    # outlet face
    Uf[:,nx+1] = U[:,nx]
# pressure correction link coefficents
def solve_pressure():
    pold = copy.deepcopy(p_prime)
#Bottom wall Block 2
    i = 1
    j = nx1+1

    ap_e = 0.5 * rho * dy * dy * (1 / ap[i][j] + 1 / ap[i][j+1])
    ap_w = 0
    ap_n = 0.5 * rho * dx * dx * (1 / ap[i][j] + 1 / ap[i+1][j])
    ap_s = 0
    ap_o = ap_e + ap_w + ap_n + ap_s
    Sp = -rho * ((Uf[i][j+1] - Uf[i][j]) * dy + (Vf[i+1][j] - Vf[i][j]) * dx)
    p_prime[i][j] = (Sp + ap_e * p_prime[i][j+1] + ap_w * p_prime[i ][j-1] + ap_n * p_prime[i+1][j] + ap_s *p_prime[i-1][j]) / ap_o
    p_prime[i][j] = p_primeold[i][j] + alpha_prime*(p_prime[i][j] - p_primeold[i][j])

    for j in range(nx1+2,nx):
        ap_e = 0.5*rho*dy*dy*(1 / ap[i][j] + 1 / ap[i][j + 1])
        ap_w = 0.5*rho*dy*dy*(1 / ap[i][j] + 1 / ap[i][j - 1])
        ap_n = 0.5*rho*dx*dx*(1 / ap[i][j] + 1 / ap[i + 1][j])
        ap_s = 0
        ap_o = ap_e + ap_w + ap_n + ap_s
        Sp = -rho*((Uf[i][j + 1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)
        p_prime[i][j] = (Sp + ap_e*p_prime[i][j+1] + ap_w*p_prime[i][j-1] + ap_n*p_prime[i+1][j] + ap_s*p_prime[i-1][j]) / ap_o
        p_prime[i][j] = p_primeold[i][j] + alpha_prime * (p_prime[i][j] - p_primeold[i][j])
    j = nx
    ap_e = 0
    ap_w = 0.5*rho*dy*dy*(1 / ap[i][j - 1])
    ap_n = 0.5*rho*dx*dx*(1 / ap[i][j] + 1 / ap[i + 1][j]);
    ap_s = 0
    ap_o = ap_e + ap_w + ap_n + ap_s + rho * dy * dy / ap[i][j]
    Sp = -rho*((Uf[i][j + 1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)

    p_prime[i][j] = (Sp + ap_e*p_prime[i][j+1] + ap_w*p_prime[i][j-1] + ap_n*p_prime[i+1][j] + ap_s*p_prime[i-1][j])/ap_o
    p_prime[i][j] = p_primeold[i][j] + alpha_prime*(p_prime[i][j] - p_primeold[i][j])

    # vertical wall
    for i in range (2 ,ny1+1):
        j = nx1 + 1
        ap_e = 0.5 * rho * dy * dy * (1 / ap[i][j] + 1 / ap[i][j + 1])
        ap_w = 0
        ap_n = 0.5 * rho * dx * dx * (1 / ap[i][j] + 1 / ap[i + 1][j])
        ap_s = 0.5 * rho * dx * dx * (1 / ap[i][j] + 1 / ap[i - 1][j])
        ap_o = ap_e + ap_w + ap_n + ap_s
        Sp = -rho * ((Uf[i][j + 1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)
        p_prime[i][j] = (Sp + ap_e * p_prime[i][j+1] + ap_w * p_prime[i][j-1] + ap_n * p_prime[i+1][j] + ap_s *p_prime[i-1][j])/ ap_o
        p_prime[i][j] = p_primeold[i][j] + alpha_prime * (p_prime[i][j] - p_primeold[i][j])
        for j in range(nx1+2,nx):
                ap_e = 0.5*rho*dy*dy*(1 / ap[i][j] + 1 / ap[i][j + 1])
                ap_w = 0.5*rho*dy*dy*(1 / ap[i][j] + 1 / ap[i][j - 1])
                ap_n = 0.5*rho*dx*dx*(1 / ap[i][j] + 1 / ap[i + 1][j])
                ap_s = 0.5*rho*dx*dx*(1 / ap[i][j] + 1 / ap[i - 1][j])
                ap_o = ap_e + ap_w + ap_n + ap_s
                Sp = -rho*((Uf[i][j + 1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)
                p_prime[i][j] = (Sp + ap_e*p_prime[i][j+1] + ap_w*p_prime[i][j-1] + ap_n*p_prime[i+1][j] + ap_s*p_prime[i-1][j]) / ap_o
                p_prime[i][j] = p_primeold[i][j] + alpha_prime * (p_prime[i][j] - p_primeold[i][j])
        j = nx

        ap_e = 0
        ap_w = 0.5 * rho * dy * dy * (1 / ap[i][j - 1])
        ap_n = 0.5 * rho * dx * dx * (1 / ap[i][j] + 1 / ap[i + 1][j])
        ap_s = 0.5 * rho * dx * dx * (1 / ap[i][j] + 1 / ap[i - 1][j])
        ap_o = ap_e + ap_w + ap_n + ap_s + rho * dy * dy / ap[i][j]
        Sp = -rho * ((Uf[i][j + 1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)
        p_prime[i][j] = (Sp + ap_e * p_prime[i][j+1] + ap_w * p_prime[i][j-1] + ap_n * p_prime[i+1][j] + ap_s *p_prime[i-1][j]) /ap_o
        p_prime[i][j] = p_primeold[i][j] + alpha_prime*(p_prime[i][j] - p_primeold[i][j])
    i = ny1 + 1
    j = 1

    ap_e = 0.5 * rho * dy * dy * (1 / ap[i][j] + 1 / ap[i][j + 1])
    ap_w = 0
    ap_n = 0.5 * rho * dx * dx * (1 / ap[i][j] + 1 / ap[i + 1][j])
    ap_s = 0
    ap_o = ap_e + ap_w + ap_n + ap_s
    Sp = -rho * ((Uf[i][j + 1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)
    p_prime[i][j] = (Sp + ap_e * p_prime[i][j+1] + ap_w * p_prime[i][j-1] + ap_n * p_prime[i+1][j] + ap_s*p_prime[i-1][j]) / ap_o
    p_prime[i][j] = p_primeold[i][j] + alpha_prime*(p_prime[i][j] - p_primeold[i][j])

    for j in range(2, nx1+1):
        ap_e = 0.5*rho*dy*dy*(1 / ap[i][j] + 1 / ap[i][j + 1])
        ap_w = 0.5*rho*dy*dy*(1 / ap[i][j] + 1 / ap[i][j - 1])
        ap_n = 0.5*rho*dx*dx*(1 / ap[i][j] + 1 / ap[i + 1][j])
        ap_s = 0
        ap_o = ap_e + ap_w + ap_n + ap_s
        Sp = -rho*((Uf[i][j + 1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)
        p_prime[i][j] = (Sp + ap_e*p_prime[i][j+1] + ap_w*p_prime[i][j-1] + ap_n*p_prime[i+1][j] + ap_s*p_prime[i-1][j]) / ap_o
        p_prime[i][j] = p_primeold[i][j] + alpha_prime * (p_prime[i][j] - p_primeold[i][j])

    for j in range(nx1+1,nx):
        ap_e = 0.5*rho*dy*dy*(1 / ap[i][j] + 1 / ap[i][j + 1])
        ap_w = 0.5*rho*dy*dy*(1 / ap[i][j] + 1 / ap[i][j - 1])
        ap_n = 0.5*rho*dx*dx*(1 / ap[i][j] + 1 / ap[i + 1][j])
        ap_s = 0.5*rho*dx*dx*(1 / ap[i][j] + 1 / ap[i - 1][j])
        ap_o = ap_e + ap_w + ap_n + ap_s
        Sp = -rho*((Uf[i][j + 1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)
        p_prime[i][j] = (Sp + ap_e*p_prime[i][j+1] + ap_w*p_prime[i][j+1] + ap_n*p_prime[i+1][j] + ap_s*p_prime[i-1][j]) / ap_o
        p_prime[i][j] = p_primeold[i][j] + alpha_prime * (p_prime[i][j] - p_primeold[i][j])
    j = nx

    ap_e = 0
    ap_w = 0.5*rho*dy*dy*(1 / ap[i][j - 1])
    ap_n = 0.5*rho*dx*dx*(1 / ap[i][j] + 1 / ap[i + 1][j])
    ap_s = 0.5*rho*dx*dx*(1 / ap[i][j] + 1 / ap[i - 1][j])
    ap_o = ap_e + ap_w + ap_n + ap_s + rho * dy * dy / ap[i][j]
    Sp = -rho*((Uf[i][j+1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)
    p_prime[i][j] = (Sp + ap_e*p_prime[i][j+1] + ap_w*p_prime[i][j-1] + ap_n*p_prime[i+1][j] + ap_s*p_prime[i-1][j]) / ap_o
    p_prime[i][j] = p_primeold[i][j] + alpha_prime*(p_prime[i][j] - p_primeold[i][j])

    for i in range(ny1+2,ny):
        j = 1

        ap_e = 0.5 * rho * dy * dy * (1 / ap[i][j] + 1 / ap[i][j + 1])
        ap_w = 0
        ap_n = 0.5 * rho * dx * dx * (1 / ap[i][j] + 1 / ap[i + 1][j])
        ap_s = 0.5 * rho * dx * dx * (1 / ap[i][j] + 1 / ap[i - 1][j])
        ap_o = ap_e + ap_w + ap_n + ap_s
        Sp = -rho * ((Uf[i][j + 1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)

        p_prime[i][j] = (Sp + ap_e * p_prime[i][j+1] + ap_w * p_prime[i][j-1] + ap_n * p_prime[i+1][j] + ap_s * p_prime[i-1][j]) / ap_o
        p_prime[i][j] = p_primeold[i][j] + alpha_prime * (p_prime[i][j] - p_primeold[i][j])
        for j in range (2, nx):
            ap_e = 0.5 * rho * dy * dy * (1 / ap[i][j] + 1 / ap[i][j + 1])
            ap_w = 0.5 * rho * dy * dy * (1 / ap[i][j] + 1 / ap[i][j - 1])
            ap_n = 0.5 * rho * dx * dx * (1 / ap[i][j] + 1 / ap[i + 1][j])
            ap_s = 0.5 * rho * dx * dx * (1 / ap[i][j] + 1 / ap[i - 1][j])
            ap_o = ap_e + ap_w + ap_n + ap_s
            Sp = -rho * ((Uf[i][j + 1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)
            p_prime[i][j] = (Sp + ap_e * p_prime[i][j+1] + ap_w * p_prime[i][j-1] + ap_n * p_prime[i+1][j] + ap_s *p_prime[i-1][j]) / ap_o
            p_prime[i][j] = p_primeold[i][j] + alpha_prime * (p_prime[i][j] - p_primeold[i][j])

        j = nx

        ap_e = 0
        ap_w = 0.5 * rho * dy * dy * (1 / ap[i][j])
        ap_n = 0.5 * rho * dx * dx * (1 / ap[i][j] + 1 / ap[i + 1][j])
        ap_s = 0.5 * rho * dx * dx * (1 / ap[i][j] + 1 / ap[i - 1][j])
        ap_o = ap_e + ap_w + ap_n + ap_s + rho * dy * dy / ap[i][j]
        Sp = -rho * ((Uf[i][j+1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)

        p_prime[i][j] = (Sp + ap_e * p_prime[i][j+1] + ap_w * p_prime[i][j-1] + ap_n * p_prime[i+1][j] + ap_s *p_prime[i-1][j]) / ap_o
        p_prime[i][j] = p_primeold[i][j] + alpha_prime * (p_prime[i][j] - p_primeold[i][j])

    i = ny
    j = 1
    ap_e = 0.5 * rho * dy * dy * (1 / ap[i][j] + 1 / ap[i][j + 1])
    ap_w = 0
    ap_n = 0
    ap_s = 0.5 * rho * dx * dx * (1 / ap[i][j] + 1 / ap[i - 1][j])
    ap_o = ap_e + ap_w + ap_n + ap_s
    Sp = -rho * ((Uf[i][j + 1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)

    p_prime[i][j] = (Sp + ap_e * p_prime[i][j+1] + ap_w * p_prime[i][j-1] + ap_n * p_prime[i+1][j] + ap_s *p_prime[i-1][j]) / ap_o
    p_prime[i][j] = p_primeold[i][j] + alpha_prime*(p_prime[i][j] - p_primeold[i][j])

    for j in range(2,nx):
        ap_e = 0.5*rho*dy*dy*(1 / ap[i][j] + 1 / ap[i][j + 1])
        ap_w = 0.5*rho*dy*dy*(1 / ap[i][j] + 1 / ap[i][j - 1])
        ap_n = 0
        ap_s = 0.5*rho*dx*dx*(1 / ap[i][j] + 1 / ap[i - 1][j])
        ap_o = ap_e + ap_w + ap_n + ap_s
        Sp = -rho*((Uf[i][j + 1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)
        p_prime[i][j] = (Sp + ap_e*p_prime[i][j+1] + ap_w*p_prime[i][j-1] + ap_n*p_prime[i+1][j] + ap_s*p_prime[i-1][j]) / ap_o
        p_prime[i][j] = p_primeold[i][j] + alpha_prime * (p_prime[i][j] - p_primeold[i][j])
    j = nx
    ap_e = 0
    ap_w = 0.5*rho*dy*dy*(1 / ap[i][j - 1])
    ap_n = 0
    ap_s = 0.5*rho*dx*dx*(1 / ap[i][j] + 1 / ap[i - 1][j])
    ap_o = ap_e + ap_w + ap_n + ap_s + rho * dy * dy / ap[i][j]
    Sp = -rho*((Uf[i][j+1] - Uf[i][j]) * dy + (Vf[i + 1][j] - Vf[i][j]) * dx)
    p_prime[i][j] = (Sp + ap_e*p_prime[i][j+1] + ap_w*p_prime[i][j-1] + ap_n*p_prime[i+1][j] + ap_s*p_prime[i-1][j]) / ap_o
    p_prime[i][j] = p_primeold[i][j] + alpha_prime*(p_prime[i][j] - p_primeold[i][j])
    res_p = np.linalg.norm(abs(pold-p_prime))
    return res_p

# correction of pressure
def pressure_correction():
    for j in range (nx1+1,nx+1):
        for i in range (1,ny1+1):
            p[i][j] += alpha_p * p_prime[i][j]

    for j in range (1, nx+1):
        for i in range (ny1+1,ny+1):
            p[i][j] += alpha_p * p_prime[i][j]

def velocity_correction():
    # correction of velocities
    j = 1
    for i in range(ny1+1 ,ny+1):
        U[i][j] += alpha_vel * (p_prime[i][j] - p_prime[i][j + 1]) * 0.5 * dy / ap[i][j]

    for j in range(2,nx):
        for i in range(ny1+1,ny+1):
            U[i][j] += alpha_vel * (p_prime[i][j - 1] - p_prime[i][j + 1]) * 0.5 * dy / ap[i][j]

    j = nx
    for i in range(1,ny+1):
        U[i][j] += alpha_vel * (p_prime[i][j - 1] - p_prime[i][j]) * 0.5 * dy / ap[i][j]
    for j in range(nx1+2,nx):
        for i in range(1,ny1+1):
            U[i][j] += alpha_vel * (p_prime[i][j - 1] - p_prime[i][j + 1]) * 0.5 * dy / ap[i][j]
    j = nx1+1;
    for i in range(1,ny1+1):
        U[i][j] += alpha_vel * (p_prime[i][j] - p_prime[i][j + 1]) * 0.5 * dy / ap[i][j]

    #v correction

    i = 1
    for j in range(nx1+1,nx+1):
        V[i][j] += alpha_vel * (p_prime[i][j] - p_prime[i + 1][j]) * 0.5 * dx / ap[i][j]

    for i in range(2,ny):
        for j in range(nx1+1, nx+1):
            V[i][j] += alpha_vel * (p_prime[i - 1][j] - p_prime[i + 1][j]) * 0.5 * dx / ap[i][j]

    i = ny
    for j in range(1,nx+1):
        V[i][j] += alpha_vel * (p_prime[i - 1][j] - p_prime[i][j]) * 0.5 * dx / ap[i][j]

    for i in range(ny1+2, ny):
        for j in range(1,nx1+1):
            V[i][j] += alpha_vel * (p_prime[i - 1][j] - p_prime[i + 1][j]) * 0.5 * dx / ap[i][j]

    i = ny1+1
    for j in range(1,nx1+1):
        V[i][j] += alpha_vel * (p_prime[i][j] - p_prime[i + 1][j]) * 0.5 * dx / ap[i][j]

# face velocity corrections
def face_vel_correction():

    for j in range(nx1+2,nx+1):
        for i in range(1,ny1+1):
            Uf[i][j] += alpha_vel * (1 / ap[i][j - 1] + 1 / ap[i][j]) * 0.5 * dy * (p_prime[i][j - 1] - p_prime[i][j])

    for j in range(2,nx+1):
        for i in range(ny1+1, ny+1):
            Uf[i][j] += alpha_vel * (1 / ap[i][j - 1] + 1 / ap[i][j]) * 0.5 * dy * (p_prime[i][j - 1] - p_prime[i][j])

    #v face velocity correction

    for i in range(ny1+2,ny+1):
        for j in range(1,nx1+1):
            Vf[i][j] += alpha_vel * (1 / ap[i - 1][j] + 1 / ap[i][j]) * 0.5 * dx * (p_prime[i - 1][j] - p_prime[i][j])

    for i in range(2,ny+1):
        for j in range(nx1+1,nx+1):
            Vf[i][j] += alpha_vel * (1 / ap[i - 1][j] + 1 / ap[i][j]) * 0.5 * dx * (p_prime[i - 1][j] - p_prime[i][j])
    print(nx)
    for i in range(1,ny+1):
        Uf[i][nx+1] = U[i][nx]
# SIMPLE

for k in range(0,20):
    link_coefficients()
    res_vel = solve_xandy_momentum()
    face_velocities()
    res_pr = solve_pressure()
    print(res_vel,res_pr)
    pressure_correction()
    velocity_correction()
    face_vel_correction()
plt.contourf(X,Y,U,100,cmap = 'RdGy')
plt.colorbar()
plt.show()

















































