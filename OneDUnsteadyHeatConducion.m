clear; clc;

%%This code is to solve 1D Unsteady Heat Conduction Problem by Finite
%%Difference Method, Implicit Scheme:

prompt='Enter the length of the bar: ';
L = input(prompt);
prompt='Enter the value of thermal diffusivity: ';
Alpha = input(prompt);
prompt='Enter the value of temperature at the first boundary: ';
T0 = input(prompt);
prompt='Enter the value of temperature at the second boundary: ';
T_end = input(prompt);
prompt='Enter the value of the initial temperature of the bar: ';
T_i = input(prompt);
prompt='Enter the number of divisions: ';
n = input(prompt);
prompt='Enter the final time: ';
t_f = input(prompt);
prompt='Enter the number of time steps: ';
n_t = input(prompt);

%%Calculating dx and dt
dx=L/n;
dt=t_f/n_t;

%%Calculating the coefficients
q=1+(2*Alpha*dt)/(dx)^2;
p=Alpha*dt/(dx)^2;

%%Filling the coeff. matrix
for i=1:n-1
    for j=1:n-1
        A(i,j)=0;
        if i==j
            A(i,j)=q;
        end
        if abs(j-i)==1
            A(i,j)=-p;
        end
    end
end

%%Giving the initial Temperature of each interior point
for k=1:n-1
    T(k,1)=T_i;
end

%%Filling the constant vector for the first iteration
b(1)=p*T0;
b(n-1,1)=p*T_end;
for u=2:n-2
    b(u)=T(u);
end

x_axis=linspace(0,L,n+1);

%%Solving Ax=b and updating the constant vector for each time step
for t=1:n_t
     x=A\b;
     T(1:n-1,t)=x;

     TT(1,t)=T0;
     TT(n+1,t)=T_end;
     TT(2:n,t)=T(1:n-1,t);

     b(1,1)=p*T0+x(1,1);
     b(n-1,1)=p*T_end+x(n-1,1);
     b(2:n-2,1)=x(2:n-2,1);

     plot (x_axis,TT(:,t)); hold on;
end