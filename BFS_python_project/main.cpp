#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>

using namespace std;

int nx, ny, nx_in, ny_in;
double dx, dy, ar, ra, length_in, length_out, height_in, height_out, length;
double rho, mu, g;

void init_ps(double** p, double** p_prime){
    int i,j;
    for(i=0; i<nx+2; i++)
    {
        for(j=0; j<ny+2; j++)
        {
            p[i][j] = p_prime[i][j] = 0;
        }
    }
}

void init_us(double** u, double** u_hat, double** u_face, double u_in){
    int i,j;
    for(i=0; i<nx+2; i++)
    {
        for(j=0; j<ny+2; j++)
        {
            u[i][j] = u_hat[i][j] = u_face[i][j] = u_in;
        }
    }
    for(j=0; j<ny+2; j++)
    {
        u_face[nx+2][j] = u_in;
    }

    i = nx_in+1;
    for(j=1; j<ny_in+1; j++)
    {
        u_face[i][j] = 0;
    }
}

void init_vs(double** v, double** v_hat, double** v_face){
    int i,j;
    for(i=0; i<nx+2; i++)
    {
        for(j=0; j<ny+2; j++)
        {
           v[i][j] = v_hat[i][j] = v_face[i][j] = 0;
        }
    }
    for(i=0; i<nx+2; i++)
    {
        v_face[i][ny+2] = 0;
    }
}

void init_ghost_cells(double** u, double** u_hat)
{
    int i,j;

    j = ny + 1;
    for(i=1; i<nx+2; i++)
    {
        u[i][j] = 0;
        u_hat[i][j] = 0;
    }

    j = 0;
    for(i=nx_in; i<nx+2; i++)
    {
        u[i][j] = 0;
        u_hat[i][j] = 0;
    }

    i = nx_in;
    for(j = 1; j<ny_in+1; j++)
    {
        u[i][j] = 0;
        u_hat[i][j] = 0;
    }

    j = ny_in;
    for(i=1; i<nx_in; i++)
    {
        u[i][j] = 0;
        u_hat[i][j] = 0;
    }

    for(i=0; i<nx_in+1; i++)
    {
        for(j=0; j<ny_in+1; j++)
        {
            u[i][j] = 0;
            u_hat[i][j] = 0;
        }
    }
}

void calc_mom_links(double** ao, double** ae, double** aw, double** an, double** as, double** Sx, double** Sy,double** u, double** v, double** u_face, double** v_face, double** p, double alpha_uv)
{
    double De, Dw, Dn, Ds, Fe, Fw, Fn, Fs;
    int i,j;

    De = Dw = mu*ar;
    Dn = Ds = mu*ra;

    // interior block 1;
    for(i=2; i<nx_in+1; i++)
    {
        for(j=ny_in+2; j<ny; j++)
        {
            Fe = rho*dy*u_face[i+1][j];
            Fw = rho*dy*u_face[i][j];
            Fn = rho*dx*v_face[i][j+1];
            Fs = rho*dx*v_face[i][j];

            ae[i][j] = De + max(0.0, -Fe);
            aw[i][j] = Dw + max(0.0,  Fw);
            an[i][j] = Dn + max(0.0, -Fn);
            as[i][j] = Ds + max(0.0,  Fs);
            ao[i][j] = De + Dw + Dn + Ds + max(0.0, Fe) + max(0.0, -Fw) + max(0.0, Fn) + max(0.0, -Fs);
            Sx[i][j] = (0.5*(p[i-1][j] - p[i+1][j])*dy);
            Sy[i][j] = 0.5*(p[i][j-1] - p[i][j+1])*dx - rho*g*dx*dy;
        }
    }

    // Second Interior block;
    for(i=nx_in+2; i<nx; i++)
    {
        for(j=2; j<ny; j++)
        {  
            Fe = rho*dy*u_face[i+1][j];
            Fw = rho*dy*u_face[i][j];
            Fn = rho*dx*v_face[i][j+1];
            Fs = rho*dx*v_face[i][j];

            ae[i][j] = De + max(0.0, -Fe);
            aw[i][j] = Dw + max(0.0,  Fw);
            an[i][j] = Dn + max(0.0, -Fn);
            as[i][j] = Ds + max(0.0,  Fs);
            ao[i][j] = De + Dw + Dn + Ds + max(0.0, Fe) + max(0.0, -Fw) + max(0.0, Fn) + max(0.0, -Fs);
            Sx[i][j] = 0.5*(p[i-1][j] - p[i+1][j])*dy;
            Sy[i][j] = 0.5*(p[i][j-1] - p[i][j+1])*dx - rho*g*dx*dy;
        }
    }

    // Left Inlet boundary
    i = 1;
    for(j=ny_in+2; j<ny; j++)
    {
        Fe = rho*dy*u_face[i+1][j];
        Fw = rho*dy*u_face[i][j];
        Fn = rho*dx*v_face[i][j+1];
        Fs = rho*dx*v_face[i][j];

        ae[i][j] = De + max(0.0, -Fe);
        aw[i][j] = 0;
        an[i][j] = Dn + max(0.0, -Fn);
        as[i][j] = Ds + max(0.0,  Fs);
        ao[i][j] = De + Dn + Ds + max(0.0, Fe) + max(0.0, Fn) + max(0.0, -Fs) + 2*Dw;
        Sx[i][j] = 0.5*(p[i][j] - p[i+1][j])*dy + rho*u_face[i][j]*u_face[i][j]*dy + 2*Dw*u_face[i][j];
        Sy[i][j] = 0.5*(p[i][j-1] - p[i][j+1])*dx - rho*g*dx*dy;
    }

    // Top  wall boundary
    j = ny;
    for(i=2; i<nx; i++)
    {
        Fe = rho*dy*u_face[i+1][j];
        Fw = rho*dy*u_face[i][j];
        Fn = rho*dx*v_face[i][j+1];
        Fs = rho*dx*v_face[i][j];

        ae[i][j] = De + max(0.0, -Fe);
        aw[i][j] = Dw + max(0.0,  Fw);
        an[i][j] = 0;
        as[i][j] = Ds + max(0.0,  Fs);
        ao[i][j] = De + Dw + Ds + max(0.0, Fe) + max(0.0, -Fw) + max(0.0, -Fs) + 2*Dn;
        Sx[i][j] = 0.5*(p[i-1][j] - p[i+1][j])*dy;
        Sy[i][j] = 0.5*(p[i][j-1] - p[i][j])*dx - rho*g*dx*dy;
    }

    // Bottom wall Block 1
    j = ny_in+1;
    for(i=2; i<nx_in+1; i++)
    {
        Fe = rho*dy*u_face[i+1][j];
        Fw = rho*dy*u_face[i][j];
        Fn = rho*dx*v_face[i][j+1];
        Fs = rho*dx*v_face[i][j];

        ae[i][j] = De + max(0.0, -Fe);
        aw[i][j] = Dw + max(0.0,  Fw);
        an[i][j] = Dn + max(0.0, -Fn);
        as[i][j] = 0;
        ao[i][j] = De + Dw + Dn + max(0.0, Fe) + max(0.0, -Fw) + max(0.0, Fn) + 2*Ds;
        Sx[i][j] = 0.5*(p[i-1][j] - p[i+1][j])*dy;
        Sy[i][j] = 0.5*(p[i][j] - p[i][j+1])*dx - rho*g*dx*dy;
    }

    // Left wall Block 2
    i = nx_in+1;
    for(j=2; j<ny_in+1; j++)
    {
        Fe = rho*dy*u_face[i+1][j];
        Fw = rho*dy*u_face[i][j];
        Fn = rho*dx*v_face[i][j+1];
        Fs = rho*dx*v_face[i][j];

        ae[i][j] = De + max(0.0, -Fe);
        aw[i][j] = 0;
        an[i][j] = Dn + max(0.0, -Fn);
        as[i][j] = Ds + max(0.0,  Fs);
        ao[i][j] = De + Dn + Ds + max(0.0, Fe) + max(0.0, Fn) + max(0.0, -Fs) + 2*Dw;
        Sx[i][j] = 0.5*(p[i][j] - p[i+1][j])*dy;
        Sy[i][j] = 0.5*(p[i][j-1] - p[i][j+1])*dx - rho*g*dx*dy;
    }

    // Bottom wall Block 2
    j = 1;
    for(i=nx_in+2; i<nx; i++)
    {
        Fe = rho*dy*u_face[i+1][j];
        Fw = rho*dy*u_face[i][j];
        Fn = rho*dx*v_face[i][j+1];
        Fs = rho*dx*v_face[i][j];

        ae[i][j] = De + max(0.0, -Fe);
        aw[i][j] = Dw + max(0.0,  Fw);
        an[i][j] = Dn + max(0.0, -Fn);
        as[i][j] = 0;
        ao[i][j] = De + Dw + Dn + max(0.0, Fe) + max(0.0, -Fw) + max(0.0, Fn) + 2*Ds;
        Sx[i][j] = 0.5*(p[i-1][j] - p[i+1][j])*dy;
        Sy[i][j] = 0.5*(p[i][j] - p[i][j+1])*dx - rho*g*dx*dy;
    }

    // Outlet
    i = nx;
    for(j=2; j<ny; j++)
    {
        Fe = rho*dy*u_face[i+1][j];
        Fw = rho*dy*u_face[i][j];
        Fn = rho*dx*v_face[i][j+1];
        Fs = rho*dx*v_face[i][j];

        ae[i][j] = 0;
        aw[i][j] = Dw + max(0.0,  Fw);
        an[i][j] = Dn + max(0.0, -Fn);
        as[i][j] = Ds + max(0.0,  Fs);
        ao[i][j] = Dw + Dn + Ds + max(0.0, Fe) + max(0.0, -Fw) + max(0.0, Fn) + max(0.0, -Fs);
        Sx[i][j] = 0.5*(p[i-1][j] + p[i][j])*dy;
        Sy[i][j] = 0.5*(p[i][j-1] - p[i][j+1])*dx - rho*g*dx*dy;
    }

    // Interior Cells b/w Blocks 1 and 2
    i = nx_in+1;
    for(j=ny_in+1; j<ny; j++)
    {
        Fe = rho*dy*u_face[i+1][j];
        Fw = rho*dy*u_face[i][j];
        Fn = rho*dx*v_face[i][j+1];
        Fs = rho*dx*v_face[i][j];

        ae[i][j] = De + max(0.0, -Fe);
        aw[i][j] = Dw + max(0.0,  Fw);
        an[i][j] = Dn + max(0.0, -Fn);
        as[i][j] = Ds + max(0.0,  Fs);
        ao[i][j] = De + Dw + Dn + Ds + max(0.0, Fe) + max(0.0, -Fw) + max(0.0, Fn) + max(0.0, -Fs);
        Sx[i][j] = 0.5*(p[i-1][j] - p[i+1][j])*dy;
        Sy[i][j] = 0.5*(p[i][j-1] - p[i][j+1])*dx - rho*g*dx*dy;
    }

    // Inlet cell at top wall
    i = 1;
    j = ny;
    
    Fe = rho*dy*u_face[i+1][j];
    Fw = rho*dy*u_face[i][j];
    Fn = rho*dx*v_face[i][j+1];
    Fs = rho*dx*v_face[i][j];

    ae[i][j] = De + max(0.0, -Fe);
    aw[i][j] = 0;
    an[i][j] = 0;
    as[i][j] = Ds + max(0.0,  Fs);
    ao[i][j] = De + Ds + max(0.0, Fe) + max(0.0, -Fs) + 2*Dn + 2*Dw;
    Sx[i][j] = 0.5*(p[i][j] - p[i+1][j])*dy + rho*u_face[i][j]*u_face[i][j]*dy + 2*Dw*u_face[i][j];
    Sy[i][j] = 0.5*(p[i][j-1] - p[i][j])*dx - rho*g*dx*dy;

    // Inlet cell at Bottom wall
    i = 1;
    j = ny_in+1;

    Fe = rho*dy*u_face[i+1][j];
    Fw = rho*dy*u_face[i][j];
    Fn = rho*dx*v_face[i][j+1];
    Fs = rho*dx*v_face[i][j];

    ae[i][j] = De + max(0.0, -Fe);
    aw[i][j] = 0;
    an[i][j] = Dn + max(0.0, -Fn);
    as[i][j] = 0;
    ao[i][j] = De + Dn + max(0.0, Fe) + max(0.0, Fn) + 2*Ds + 2*Dw;
    Sx[i][j] = 0.5*(p[i][j] - p[i+1][j])*dy + rho*u_face[1][j]*u_face[1][j]*dy + 2*Dw*u_face[1][j];
    Sy[i][j] = 0.5*(p[i][j] - p[i][j+1])*dx - rho*g*dx*dy;

    // Outlet cell at Top wall
    i = nx;
    j = ny;
    
    Fe = rho*dy*u_face[i+1][j];
    Fw = rho*dy*u_face[i][j];
    Fn = rho*dx*v_face[i][j+1];
    Fs = rho*dx*v_face[i][j];

    ae[i][j] = 0;
    aw[i][j] = Dw + max(0.0,  Fw);
    an[i][j] = 0;
    as[i][j] = Ds + max(0.0,  Fs);
    ao[i][j] = Dw + Ds + max(0.0, Fe) + max(0.0, -Fw) + max(0.0, -Fs) + 2*Dn;
    Sx[i][j] = 0.5*(p[i-1][j] + p[i][j])*dy;
    Sy[i][j] = 0.5*(p[i][j-1] - p[i][j])*dx - rho*g*dx*dy;

    // Outlet cell at Bottom wall
    i = nx;
    j = 1;
    
    Fe = rho*dy*u_face[i+1][j];
    Fw = rho*dy*u_face[i][j];
    Fn = rho*dx*v_face[i][j+1];
    Fs = rho*dx*v_face[i][j];

    ae[i][j] = 0;
    aw[i][j] = Dw + max(0.0,  Fw);
    an[i][j] = Dn + max(0.0, -Fn);
    as[i][j] = 0;
    ao[i][j] = Dw + Dn + max(0.0, Fe) + max(0.0, -Fw) + max(0.0, Fn) + 2*Ds;
    Sx[i][j] = 0.5*(p[i-1][j] + p[i][j])*dy;
    Sy[i][j] = 0.5*(p[i][j] - p[i][j+1])*dx - rho*g*dx*dy;
   
   // Bottom left cell Block 2
   i = nx_in+1;
   j = 1;

    Fe = rho*dy*u_face[i+1][j];
    Fw = rho*dy*u_face[i][j];
    Fn = rho*dx*v_face[i][j+1];
    Fs = rho*dx*v_face[i][j];

    ae[i][j] = De + max(0.0, -Fe);
    aw[i][j] = 0;
    an[i][j] = Dn + max(0.0, -Fn);
    as[i][j] = 0;
    ao[i][j] = De + Dn + max(0.0, Fe) + max(0.0, Fn) + 2*Ds + 2*Dw;
    Sx[i][j] = 0.5*(p[i][j] - p[i+1][j])*dy;
    Sy[i][j] = 0.5*(p[i][j] - p[i][j+1])*dx - rho*g*dx*dy;

    for(i=1; i<nx+1; i++)
    {
        for(j=1; j<ny+1; j++)
        {
            ae[i][j] = alpha_uv*ae[i][j];
            aw[i][j] = alpha_uv*aw[i][j];
            an[i][j] = alpha_uv*an[i][j];
            as[i][j] = alpha_uv*as[i][j];
            Sx[i][j] = alpha_uv*Sx[i][j] + (1-alpha_uv)*ao[i][j]*u[i][j];
            Sy[i][j] = alpha_uv*Sy[i][j] + (1-alpha_uv)*ao[i][j]*v[i][j];
        }
    }
}

void solve_velocity(double** ao, double** ae, double** aw, double** an, double** as, double** S,double** vel_hat, double** vel_hat_old, double**vel,  double alpha)
{
    int i,j,k=0;
    int gs_iter = 200;
    double norm = 1, error;
    double residue = 0;

    for(i=1; i<nx+1; i++)
    {
        for(j=1; j<ny+1; j++)
        {
            residue += fabs(S[i][j] + ae[i][j]*vel[i+1][j] + aw[i][j]*vel[i-1][j] + an[i][j]*vel[i][j+1] + as[i][j]*vel[i][j-1] - ao[i][j]*vel[i][j]);
        }
    }

    for(k=0; k<gs_iter; k++)
    {
        norm = 0;
        k++;
        for(i=1; i<nx+1; i++)
        {
            for(j=1; j<ny+1; j++)
            {
                vel_hat_old[i][j] = vel_hat[i][j];
            }
        }

        for(j=1; j<ny_in+1; j++)
        {
            for(i=nx_in+1; i<nx+1; i++)
            {
                vel_hat[i][j] = alpha*(S[i][j] + ae[i][j]*vel_hat[i+1][j] + aw[i][j]*vel_hat[i-1][j] + an[i][j]*vel_hat[i][j+1] + as[i][j]*vel_hat[i][j-1])/ao[i][j] + (1-alpha)*vel_hat[i][j];
            }
        }
        
        for(j=ny_in+1; j<ny+1; j++)
        {
            for(i=nx; i>0; i--)
            {
                vel_hat[i][j] = alpha*(S[i][j] + ae[i][j]*vel_hat[i+1][j] + aw[i][j]*vel_hat[i-1][j] + an[i][j]*vel_hat[i][j+1] + as[i][j]*vel_hat[i][j-1])/ao[i][j] + (1-alpha)*vel_hat[i][j];
            }
        }

        for(i=1; i<nx+1; i++)
        {
            for(j=1; j<ny+1; j++)
            {
                error = fabs(vel_hat[i][j] - vel_hat_old[i][j]);
                if(error > norm)
                {
                    norm = error;
                }
            }
        }

    }
    cout << "uv_norm: " << norm << endl;
    // cout << "Residue: " << residue << endl;
}

void face_velocity(double** u_face, double** v_face, double** u_hat, double** v_hat, double** p, double** ao, double alpha_uv)
{
    int i,j;

    // u face velocities
    //
    // internal faces
    for(i=nx_in+3; i<nx; i++)
    {
        for(j=1; j<ny_in+1; j++)
        {
            u_face[i][j] = 0.5*(u_hat[i][j]+u_hat[i-1][j]) + 0.25*alpha_uv*dy*((p[i][j]-p[i-2][j])/ao[i-1][j] + (p[i+1][j]-p[i-1][j])/ao[i][j]) - 0.5*alpha_uv*dy*(1/ao[i][j]+1/ao[i-1][j])*(p[i][j]-p[i-1][j]);
        }
    }

    for(i=3; i<nx; i++)
    {
        for(j=ny_in+1; j<ny+1; j++)
        {
            u_face[i][j] = 0.5*(u_hat[i][j]+u_hat[i-1][j]) + 0.25*alpha_uv*dy*((p[i][j]-p[i-2][j])/ao[i-1][j] + (p[i+1][j]-p[i-1][j])/ao[i][j]) - 0.5*alpha_uv*dy*(1/ao[i][j]+1/ao[i-1][j])*(p[i][j]-p[i-1][j]);
        }
    }

    // adjacent to inlet boundary
    i = 2;
    for(j=ny_in+1; j<ny+1; j++)
    {
        u_face[i][j] = 0.5*(u_hat[i][j]+u_hat[i-1][j]) + 0.25*alpha_uv*dy*((p[i][j]-p[i-1][j])/ao[i-1][j] + (p[i+1][j]-p[i-1][j])/ao[i][j]) - 0.5*alpha_uv*dy*(1/ao[i][j]+1/ao[i-1][j])*(p[i][j]-p[i-1][j]);
    }

    // adjacent to left bottom wall
    i = nx_in + 2;
    for(j=1; j<ny_in+1; j++)
    {
        u_face[i][j] = 0.5*(u_hat[i][j]+u_hat[i-1][j]) + 0.25*alpha_uv*dy*((p[i][j]-p[i-1][j])/ao[i-1][j] + (p[i+1][j]-p[i-1][j])/ao[i][j]) - 0.5*alpha_uv*dy*(1/ao[i][j]+1/ao[i-1][j])*(p[i][j]-p[i-1][j]);
    }

    // adjacent to outlet
    i = nx;
    for(j=1; j<ny+1; j++)
    {
        u_face[i][j] = 0.5*(u_hat[i][j]+u_hat[i-1][j]) + 0.25*alpha_uv*dy*((p[i][j]-p[i-2][j])/ao[i-1][j] - (p[i-1][j])/ao[i][j]) - 0.5*alpha_uv*dy*(1/ao[i][j]+1/ao[i-1][j])*(p[i][j]-p[i-1][j]);
    }

    // outlet face
    i = nx+1;
    for(j=1; j<ny+1; j++)
    {
        u_face[i][j] = u_hat[i-1][j];
    }

    // v face velocities
    //
    // internal faces
    for(i=1; i<nx_in+1; i++)
    {
        for(j=ny_in+3; j<ny; j++)
        {
            v_face[i][j] = 0.5*(v_hat[i][j]+v_hat[i][j-1]) + 0.25*alpha_uv*dx*((p[i][j]-p[i][j-2])/ao[i][j-1] + (p[i][j+1]-p[i][j-1])/ao[i][j]) - 0.5*alpha_uv*dx*(1/ao[i][j]+1/ao[i][j-1])*(p[i][j]-p[i][j-1]);
        }
    }

    for(i=nx_in+1; i<nx+1; i++)
    {
        for(j=3; j<ny; j++)
        {
            v_face[i][j] = 0.5*(v_hat[i][j]+v_hat[i][j-1]) + 0.25*alpha_uv*dx*((p[i][j]-p[i][j-2])/ao[i][j-1] + (p[i][j+1]-p[i][j-1])/ao[i][j]) - 0.5*alpha_uv*dx*(1/ao[i][j]+1/ao[i][j-1])*(p[i][j]-p[i][j-1]);
        }
    }

    // adjacent to Top boundary
    j = ny;
    for(i=1; i<nx+1; i++)
    {
        v_face[i][j] = 0.5*(v_hat[i][j]+v_hat[i][j-1]) + 0.25*alpha_uv*dx*((p[i][j]-p[i][j-2])/ao[i][j-1] + (p[i][j]-p[i][j-1])/ao[i][j]) - 0.5*alpha_uv*dx*(1/ao[i][j]+1/ao[i][j-1])*(p[i][j]-p[i][j-1]);
    }

    // adjacent to bottom wall block 1
    j = ny_in + 2;
    for(i=1; i<nx_in+1; i++)
    {
        v_face[i][j] = 0.5*(v_hat[i][j]+v_hat[i][j-1]) + 0.25*alpha_uv*dx*((p[i][j]-p[i][j-1])/ao[i][j-1] + (p[i][j+1]-p[i][j-1])/ao[i][j]) - 0.5*alpha_uv*dx*(1/ao[i][j]+1/ao[i][j-1])*(p[i][j]-p[i][j-1]);
    }

    // adjacent to Bottom wall block 2
    j = 2;
    for(i=nx_in+1; i<nx+1; i++)
    {
        v_face[i][j] = 0.5*(v_hat[i][j]+v_hat[i][j-1]) + 0.25*alpha_uv*dx*((p[i][j]-p[i][j-1])/ao[i][j-1] + (p[i][j+1]-p[i][j-1])/ao[i][j]) - 0.5*alpha_uv*dx*(1/ao[i][j]+1/ao[i][j-1])*(p[i][j]-p[i][j-1]);
    }
}

void calc_press_links(double** ap_o, double** ap_e, double** ap_w, double** ap_n, double** ap_s, double** Sp, double** ao, double** u_face, double** v_face, double alpha_uv)
{
    int i,j;

    // Interior 1
    for(i=2; i<nx_in+1; i++)
    {
        for(j=ny_in+2; j<ny; j++)
        {
            ap_e[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i+1][j]);
            ap_w[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i-1][j]);
            ap_n[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j+1]);
            ap_s[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j-1]);
            ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] + ap_s[i][j];
            Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);
        }
    }

    // Interior 2
    for(i=nx_in+2; i<nx; i++)
    {
        for(j=2; j<ny; j++)
        {
            ap_e[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i+1][j]);
            ap_w[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i-1][j]);
            ap_n[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j+1]);
            ap_s[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j-1]);
            ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] + ap_s[i][j];
            Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);
        }
    }

    // Inlet
    i = 1;
    for(j=ny_in+2; j<ny; j++)
    {
        ap_e[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i+1][j]);
        ap_w[i][j] = 0;
        ap_n[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j+1]);
        ap_s[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j-1]); 
        ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] + ap_s[i][j];
        Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);
    }

    // Top wall
    j = ny;
    for(i=2; i<nx; i++)
    {
        ap_e[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i+1][j]);
        ap_w[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i-1][j]);
        ap_n[i][j] = 0;
        ap_s[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j-1]);
        ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] + ap_s[i][j];
        Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);
    }

    // Bottom wall Block 1
    j = ny_in + 1;
    for(i=2; i<nx_in+1; i++)
    {
        ap_e[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i+1][j]);
        ap_w[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i-1][j]);
        ap_n[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j+1]);
        ap_s[i][j] = 0;
        ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] + ap_s[i][j];
        Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);
    }

    // Left wall Block 2
    i = nx_in + 1;
    for(j=2; j<ny_in+1; j++)
    {
        ap_e[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i+1][j]);
        ap_w[i][j] = 0;
        ap_n[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j+1]);
        ap_s[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j-1]);
        ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] + ap_s[i][j];
        Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);
    }

    // Bottom wall Block 2
    j = 1;
    for(i=nx_in+2; i<nx; i++)
    {
        ap_e[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i+1][j]);
        ap_w[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i-1][j]);
        ap_n[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j+1]);
        ap_s[i][j] = 0;
        ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] + ap_s[i][j];
        Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);
    }

    // Outlet
    i = nx; 
    for(j=2; j<ny; j++)
    {
        ap_e[i][j] = 0;
        ap_w[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i-1][j]);
        ap_n[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j+1]);
        ap_s[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j-1]);
        ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] + ap_s[i][j] + alpha_uv*rho*dy*dy/ao[i][j];
        Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);
    }

    // Interior Cells b/w Blocks 1 and 2
    i = nx_in+1;
    for(j=ny_in+1; j<ny; j++)
    { 
        ap_e[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i+1][j]);
        ap_w[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i-1][j]);
        ap_n[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j+1]);
        ap_s[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j-1]);
        ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] + ap_s[i][j];
        Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);
    }

    // Top wall Inlet
    j = ny;
    i = 1;

    ap_e[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i+1][j]);
    ap_w[i][j] = 0;
    ap_n[i][j] = 0;
    ap_s[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j-1]); 
    ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] + ap_s[i][j];
    Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);

    // Bottom wall Inlet
    j = ny_in + 1;
    i = 1;
    
    ap_e[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i+1][j]);
    ap_w[i][j] = 0;
    ap_n[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j+1]);
    ap_s[i][j] = 0; 
    ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] + ap_s[i][j];
    Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);

    // Top wall Outlet
    j = ny;
    i = nx;

    ap_e[i][j] = 0;
    ap_w[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i-1][j]);
    ap_n[i][j] = 0;
    ap_s[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j-1]);
    ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] +ap_s[i][j] + alpha_uv*rho*dy*dy/ao[i][j];
    Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);

    // Bottom wall outlet
    j = 1;
    i = nx;

    ap_e[i][j] = 0;
    ap_w[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i-1][j]);
    ap_n[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j+1]);
    ap_s[i][j] = 0;
    ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] +ap_s[i][j] + alpha_uv*rho*dy*dy/ao[i][j];
    Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);

    // Bottom wall left corner 
    j = 1;
    i = nx_in+1;
    
    ap_e[i][j] = 0.5*alpha_uv*rho*dy*dy*(1/ao[i][j] + 1/ao[i+1][j]);
    ap_w[i][j] = 0;
    ap_n[i][j] = 0.5*alpha_uv*rho*dx*dx*(1/ao[i][j] + 1/ao[i][j+1]);
    ap_s[i][j] = 0;
    ap_o[i][j] = ap_e[i][j] + ap_w[i][j] + ap_n[i][j] + ap_s[i][j];
    Sp[i][j] = -rho*((u_face[i+1][j] - u_face[i][j])*dy + (v_face[i][j+1] - v_face[i][j])*dx);    
}

void solve_p_prime(double** ao, double** ae, double** aw, double** an, double** as, double** S,double** p_prime, double** p_prime_old,  double alpha)
{
    int i,j,k=0;
    int gs_iter = 200000;
    double norm = 1, error;
    double residue = 0;
    alpha = 1;

    for(i=1; i<nx+1; i++)
    {
        for(j=1; j<ny+1; j++)
        {
            residue += pow(S[i][j],2);
        }
    }
    residue = sqrt(residue);

    // for(k=0; k<gs_iter; k++)
    while(norm > 1e-9)
    {
        norm = 0;
        k++;
        for(i=1; i<nx+1; i++)
        {
            for(j=1; j<ny+1; j++)
            {
                p_prime_old[i][j] = p_prime[i][j];
            }
        }

        for(j=1; j<ny_in+1; j++)
        {
            for(i=nx_in+1; i<nx+1; i++)
            {
                p_prime[i][j] = alpha*(S[i][j] + ae[i][j]*p_prime[i+1][j] + aw[i][j]*p_prime[i-1][j] + an[i][j]*p_prime[i][j+1] + as[i][j]*p_prime[i][j-1])/ao[i][j] + (1-alpha)*p_prime[i][j];
            }
        }
        
        for(j=ny_in+1; j<ny+1; j++)
        {
            for(i=nx; i>0; i--)
            {
                p_prime[i][j] = alpha*(S[i][j] + ae[i][j]*p_prime[i+1][j] + aw[i][j]*p_prime[i-1][j] + an[i][j]*p_prime[i][j+1] + as[i][j]*p_prime[i][j-1])/ao[i][j] + (1-alpha)*p_prime[i][j];
            }
        }

        for(i=1; i<nx+1; i++)
        {
            for(j=1; j<ny+1; j++)
            {
                error = fabs(p_prime[i][j] - p_prime_old[i][j]);
                if(error > norm)
                {
                    norm = error;
                }
            }
        }
    }
    cout << "p_prime norm: " << norm << endl;
    cout << "Residue: " << residue << endl;
}

void correct_p(double** p, double** p_prime, double alpha_p)
{
    int i,j;

    for(i=1; i<nx+1; i++)
    {
        for(j=1; j<ny+1; j++)
        {
            p[i][j] += alpha_p*p_prime[i][j];
        }
    }
}

void correct_uv(double** u, double** u_hat, double** u_face, double** v, double** v_hat, double** v_face, double** p_prime, double** ao, double alpha_uv)
{
    int i,j;

    // u correction

    i = 1;
    for(j=ny_in+1; j<ny+1; j++)
    {
        u[i][j] = u_hat[i][j] + alpha_uv*(p_prime[i][j] - p_prime[i+1][j])*0.5*dy/ao[i][j];
    }

    for(i=2; i<nx; i++)
    {
        for(j=ny_in+1; j<ny+1; j++)
        {
            u[i][j] = u_hat[i][j] + alpha_uv*(p_prime[i-1][j] - p_prime[i+1][j])*0.5*dy/ao[i][j];
        }
    }

    i = nx;
    for(j=1; j<ny+1; j++)
    {
        u[i][j] = u_hat[i][j] + alpha_uv*(p_prime[i-1][j] - p_prime[i][j])*0.5*dy/ao[i][j];
    }

    for(i=nx_in+2; i<nx; i++)
    {
        for(j=1; j<ny_in+1; j++)
        {
            u[i][j] = u_hat[i][j] + alpha_uv*(p_prime[i-1][j] - p_prime[i+1][j])*0.5*dy/ao[i][j];
        }
    }

    i = nx_in+1;
    for(j=1; j<ny_in+1; j++)
    {
        u[i][j] = u_hat[i][j] + alpha_uv*(p_prime[i][j] - p_prime[i+1][j])*0.5*dy/ao[i][j];
    }

    // v correction

    j = 1;
    for(i=nx_in+1; i<nx+1; i++)
    {
        v[i][j] = v_hat[i][j] + alpha_uv*(p_prime[i][j] - p_prime[i][j+1])*0.5*dx/ao[i][j];
    }

    for(j=2; j<ny; j++)
    {
        for(i=nx_in+1; i<nx+1; i++)
        {
            v[i][j] = v_hat[i][j] + alpha_uv*(p_prime[i][j-1] - p_prime[i][j+1])*0.5*dx/ao[i][j];
        }
    }

    j = ny;
    for(i=1; i<nx+1; i++)
    {
        v[i][j] = v_hat[i][j] + alpha_uv*(p_prime[i][j-1] - p_prime[i][j])*0.5*dx/ao[i][j];
    }

    for(j=ny_in+2; j<ny; j++)
    {
        for(i=1; i<nx_in+1; i++)
        {
            v[i][j] = v_hat[i][j] + alpha_uv*(p_prime[i][j-1] - p_prime[i][j+1])*0.5*dx/ao[i][j];
        }
    }

    j = ny_in+1;
    for(i=1; i<nx_in+1; i++)
    {
        v[i][j] = v_hat[i][j] + alpha_uv*(p_prime[i][j] - p_prime[i][j+1])*0.5*dx/ao[i][j];
    }

    // u_face correction

    for(i=nx_in+2; i<nx+1; i++)
    {
        for(j=1; j<ny_in+1; j++)
        {
            u_face[i][j] += alpha_uv*(1/ao[i-1][j] + 1/ao[i][j])*0.5*dy*(p_prime[i-1][j] - p_prime[i][j]);
        }
    }

    for(i=2; i<nx+1; i++)
    {
        for(j=ny_in+1; j<ny+1; j++)
        {
            u_face[i][j] += alpha_uv*(1/ao[i-1][j] + 1/ao[i][j])*0.5*dy*(p_prime[i-1][j] - p_prime[i][j]);
        }
    }

    i = nx+1;
    for(j=1; j<ny+1; j++)
    {
        u_face[i][j] = u[i-1][j];
    }

    // v_face correction

    for(j=ny_in+2; j<ny+1; j++)
    {
        for(i=1; i<nx_in+1; i++)
        {
            v_face[i][j] += alpha_uv*(1/ao[i][j-1] + 1/ao[i][j])*0.5*dx*(p_prime[i][j-1] - p_prime[i][j]);
        }
    }

    for(j=2; j<ny+1; j++)
    {
        for(i=nx_in+1; i<nx+1; i++)
        {
            v_face[i][j] += alpha_uv*(1/ao[i][j-1] + 1/ao[i][j])*0.5*dx*(p_prime[i][j-1] - p_prime[i][j]);
        }
    }

}

int main(){

    int i,j;
    double **p, **p_prime;
    double **u, **u_hat;
    double **v, **v_hat;
    double **u_face, **v_face;
    double **ao, **ae, **aw, **an, **as, **ap_o, **ap_e, **ap_w, **ap_n, **ap_s;
    double **Sx, **Sy, **Sp;
    double **p_prime_old, **u_hat_old, **v_hat_old;
    double alpha_uv = 0.7, alpha_p = 0.5;
    string fname = "output.dat";
    ofstream outfile;

    length_in  = 0.01;
    height_in  = 0.01;
    length_out = 0.1;
    height_out = 0.02;
    length = length_in + length_out;

    ny_in = 10;
    nx_in = 20;

    dy = height_in/ny_in;
    dx = length_in/nx_in;

    nx = length/dx;
    ny = height_out/dy;

    ar = dy/dx;
    ra = 1/ar;

    rho = 1.2;
    mu = 1.8e-5;
    g = 0;

    p           = new double*[nx+2];
    p_prime     = new double*[nx+2];
    p_prime_old = new double*[nx+2];
    u           = new double*[nx+2];
    u_hat       = new double*[nx+2];
    u_hat_old   = new double*[nx+2];
    v           = new double*[nx+2];
    v_hat       = new double*[nx+2];
    v_hat_old   = new double*[nx+2];
    u_face      = new double*[nx+3];
    v_face      = new double*[nx+2];
    ao          = new double*[nx+2];
    ae          = new double*[nx+2];
    aw          = new double*[nx+2];
    an          = new double*[nx+2];
    as          = new double*[nx+2];
    Sx          = new double*[nx+2];
    Sy          = new double*[nx+2];
    ap_o        = new double*[nx+2];
    ap_e        = new double*[nx+2];
    ap_w        = new double*[nx+2];
    ap_n        = new double*[nx+2];
    ap_s        = new double*[nx+2];
    Sp          = new double*[nx+2];

    for(i=0; i<nx+2; i++)
    {
        p[i]            = new double[ny+2];
        p_prime[i]      = new double[ny+2];
        p_prime_old[i]  = new double[ny+2];
        u[i]            = new double[ny+2];
        u_hat[i]        = new double[ny+2];
        u_hat_old[i]    = new double[ny+2];
        v[i]            = new double[ny+2];
        v_hat[i]        = new double[ny+2];
        v_hat_old[i]    = new double[ny+2];
        u_face[i]       = new double[ny+2];
        v_face[i]       = new double[ny+3];
        ao[i]           = new double[ny+2];
        ae[i]           = new double[ny+2];
        aw[i]           = new double[ny+2];
        an[i]           = new double[ny+2];
        as[i]           = new double[ny+2];
        Sx[i]           = new double[ny+2];
        Sy[i]           = new double[ny+2];
        ap_o[i]         = new double[ny+2];
        ap_e[i]         = new double[ny+2];
        ap_w[i]         = new double[ny+2];
        ap_n[i]         = new double[ny+2];
        ap_s[i]         = new double[ny+2];
        Sp[i]           = new double[ny+2];
    }
    u_face[nx+2] = new double[ny+2];

    double u_in = 0.015;

    init_ps(p, p_prime);
    init_us(u, u_hat, u_face, u_in); 
    init_vs(v, v_hat, v_face);
    init_ghost_cells(u, u_hat);

    int max_iters = 200;
    double tol_out = 1e-6;

    for(int k=0; k<max_iters; k++)
    {
        calc_mom_links(ao, ae, aw, an, as, Sx, Sy, u, v, u_face, v_face, p, alpha_uv);
        solve_velocity(ao, ae, aw, an, as, Sx, u_hat, u_hat_old, u, alpha_uv);
        solve_velocity(ao, ae, aw, an, as, Sy, v_hat, v_hat_old, v, alpha_uv);
        face_velocity(u_face, v_face, u_hat, v_hat, p, ao, alpha_uv);
        calc_press_links(ap_o, ap_e, ap_w, ap_n, ap_s, Sp, ao, u_face, v_face, alpha_uv);
        solve_p_prime(ap_o, ap_e, ap_w, ap_n, ap_s, Sp, p_prime, p_prime_old, alpha_p);
        correct_p(p, p_prime, alpha_p);
        correct_uv(u, u_hat, u_face, v, v_hat, v_face, p_prime, ao, alpha_uv);        
    }

    // outfile << endl << "U_hat velocity" << endl;
    // // outfile << endl << "ap" << endl;
    // for(j=ny; j>0; j--)
    // {
    //     for(i=1; i<nx+1; i++)
    //     {
    //         outfile << u_hat[i][j] << " ";
    //     }
    //     outfile << endl;
    // }

    // outfile << endl << "V_hat velocity" << endl;
    // // outfile << endl << "ae" << endl;
    // for(j=ny; j>0; j--)
    // {
    //     for(i=1; i<nx+1; i++)
    //     {
    //         outfile << v_hat[i][j] << "\t";
    //     }
    //     outfile << endl;
    // }

    // outfile << endl << "U_face velocity" << endl;
    // // outfile << endl << "as" << endl;
    // for(j=ny; j>0; j--)
    // {
    //     for(i=1; i<nx+2; i++)
    //     {
    //         outfile << u_face[i][j] << "\t";
    //     }
    //     outfile << endl;
    // }

    // outfile << endl << "V_face velocity" << endl;
    // // outfile << endl << "Sx" << endl;
    // for(j=ny+1; j>0; j--)
    // {
    //     for(i=1; i<nx+1; i++)
    //     {
    //         outfile << v_face[i][j] << "\t";
    //     }
    //     outfile << endl;
    // }

    outfile.open("u.dat");
    // outfile << endl << "aw" << endl;
    for(j=ny; j>0; j--)
    {
        for(i=1; i<nx+1; i++)
        {
            outfile << u[i][j] << " ";
        }
        outfile << endl;
    }
    outfile.close();

    outfile.open("v.dat");
    // outfile << endl << "an" << endl;
    for(j=ny; j>0; j--)
    {
        for(i=1; i<nx+1; i++)
        {
            outfile << v[i][j] << "\t";
        }
        outfile << endl;
    }
    outfile.close();

    // outfile << endl << "p_prime" << endl;
    // // outfile << endl << "Sy" << endl;
    // for(j=ny; j>0; j--)
    // {
    //     for(i=1; i<nx+1; i++)
    //     {
    //         outfile << p_prime[i][j] << "\t";
    //     }
    //     outfile << endl;
    // }

    outfile.open("p.dat");
    for(j=ny; j>0; j--)
    {
        for(i=1; i<nx+1; i++)
        {
            outfile << p[i][j] << "\t";
        }
        outfile << endl;
    }
    outfile.close();

    // outfile << endl << "Sx" << endl;
    // // outfile << endl << "ap" << endl;
    // for(j=ny; j>0; j--)
    // {
    //     for(i=1; i<nx+1; i++)
    //     {
    //         outfile << Sx[i][j] << "\t";
    //     }
    //     outfile << endl;
    // }

    // outfile << endl << "Sy" << endl;
    // // outfile << endl << "ap" << endl;
    // for(j=ny; j>0; j--)
    // {
    //     for(i=1; i<nx+1; i++)
    //     {
    //         outfile << Sy[i][j] << "\t";
    //     }
    //     outfile << endl;
    // }

    // outfile << endl << "ao" << endl;
    // for(j=ny; j>0; j--)
    // {
    //     for(i=1; i<nx+1; i++)
    //     {
    //         outfile << ao[i][j] << "\t";
    //     }
    //     outfile << endl;
    // }

    // outfile << endl << "Sp" << endl;
    // for(j=ny; j>0; j--)
    // {
    //     for(i=1; i<nx+1; i++)
    //     {
    //         outfile << Sp[i][j] << "\t";
    //     }
    //     outfile << endl;
    // }

    // outfile << endl << "ap" << endl;
    // for(j=ny; j>0; j--)
    // {
    //     for(i=1; i<nx+1; i++)
    //     {
    //         outfile << ap_o[i][j] << "\t";
    //     }
    //     outfile << endl;
    // }

    // outfile << endl << "ae" << endl;
    // for(j=ny; j>0; j--)
    // {
    //     for(i=1; i<nx+1; i++)
    //     {
    //         outfile << ap_e[i][j] << "\t";
    //     }
    //     outfile << endl;
    // }

    // outfile << endl << "aw" << endl;
    // for(j=ny; j>0; j--)
    // {
    //     for(i=1; i<nx+1; i++)
    //     {
    //         outfile << ap_w[i][j] << "\t";
    //     }
    //     outfile << endl;
    // }

    // outfile << endl << "an" << endl;
    // for(j=ny; j>0; j--)
    // {
    //     for(i=1; i<nx+1; i++)
    //     {
    //         outfile << ap_n[i][j] << "\t";
    //     }
    //     outfile << endl;
    // }

    // outfile << endl << "as" << endl;
    // for(j=ny; j>0; j--)
    // {
    //     for(i=1; i<nx+1; i++)
    //     {
    //         outfile << ap_s[i][j] << "\t";
    //     }
    //     outfile << endl;
    // }

    outfile.close();

    return 0;
}
