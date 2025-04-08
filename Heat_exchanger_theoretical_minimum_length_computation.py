import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

T_hot_in = 90.0 #°C
T_hot_out = 45 #°C
T_cold_in = 20.0 #°C
T_cold_out = 0 # ?
m_hot = 0.8 #kg/s
m_cold = 1.2 #kg/s
d_hot =  965.0 #kg/m^3 # density of the fluid
d_cold = 998.0 #kg/m^3
mu_hot = 0.00037 #kg/m s dynamic viscosity
mu_cold = 0.001002 #kg/m s
k_hot = 0.677 # W/m K thermic conductivity
k_cold = 0.606 # W/m K
D_in_in = 0.02 # m  internal diameter internal pipe
D_in_out = 0.025 # m external diameter internal pipe
D_out_in = 0.03 # m internal diameter external pipe
Cp = 4186 # J/kg K  specific heat of the water
k_tube = 390.0 # J/m K  thermic conductivity of the copper


# evaluation of the exit temperature of cold water

T_cold_out = sp.symbols('T_cold_out')
equazione = sp.Eq(m_hot*(T_hot_out-T_hot_in),-m_cold*(T_cold_out-T_cold_in))
soluzione = sp.solve(equazione,T_cold_out)
T_cold_out = soluzione[0].evalf()

A_hot = np.pi*((D_in_in/2)**2)
v_hot = m_hot/(d_hot*A_hot)

A_cold = np.pi*(((D_out_in**2)-(D_in_out**2))/4)
v_cold = m_cold/(d_cold*A_cold)

Re_hot = (d_hot*v_hot*D_in_in)/mu_hot
Re_cold = (d_cold*v_cold*D_in_out)/mu_cold

Pr_hot = (mu_hot*Cp)/k_hot
Pr_cold = (mu_cold*Cp)/k_cold

Nu_hot = 0.023*(Re_hot**0.8)*(Pr_hot**0.3)
Nu_cold = 0.023*(Re_cold**0.8)*(Pr_cold**0.3)

h_hot = (Nu_hot*k_hot)/D_in_in
h_cold = (Nu_cold*k_cold)/D_in_out

Rf = (np.log(D_in_out/D_in_in))/(2*np.pi*k_tube)

L = 3

for i in range(10) :

    L_hp = L

    As = np.pi * D_in_out * L_hp

    # evaluation of the total exchange coefficient

    h_global = sp.symbols('h_global')
    equazione_2 = sp.Eq(1 / h_global, (1 / h_hot) + (1 / h_cold) + (Rf / As))
    soluzione_2 = sp.solve(equazione_2, h_global)
    h_global_float = soluzione_2[0].evalf()
    h_global_int = int(h_global_float)

    # realization an iterative algorithm for the determination of the theoretical length

    N = 1000
    T_hot = np.full(N, T_hot_in)
    T_cold = np.full(N, T_cold_in)
    L = 0
    T_hot_out_iter = 0

    while 1 > 0:

        for _ in range(1000):
            dx = L / N
            T_cold_new = np.copy(T_cold)
            for i in range(1, N):
                Q = h_global_int * (T_hot[i - 1] - T_cold[i - 1]) * dx
                dT_cold = Q / (m_cold * Cp)
                T_cold_new[i] = T_cold[i - 1] + dT_cold

            T_hot_new = np.copy(T_hot)
            for i in range(N - 2, -1, -1):
                Q = h_global_int * (T_hot[i + 1] - T_cold_new[i + 1]) * dx
                dT_hot = -Q / (m_hot * Cp)
                T_hot_new[i] = T_hot[i + 1] + dT_hot

            T_cold, T_hot = T_cold_new, T_hot_new
            T_hot_out_iter = T_hot


        if T_hot_out_iter[0] <= 45:
            break

        L += 0.1

T_hot_out = T_hot_out_iter[0]
print('The final temperature is '+ str(T_hot_out))
print('The theoretical length of the heat exchanger is '+ str(L))

# plot dei risultati

x = np.linspace(0,L,N)
plt.plot(x, T_hot, label='Hot Fluid',color='red')
plt.plot(x, T_cold, label='Cold Fluid',color='blue')
plt.xlabel('exchanger length (m)')
plt.ylabel('Temperature (°C)')
plt.legend()
plt.grid()
plt.show()