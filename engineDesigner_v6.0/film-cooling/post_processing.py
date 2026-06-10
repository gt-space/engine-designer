# post_processing.py
import matplotlib.pyplot as plt

def plot_results(data_log):
    # Extract arrays
    x_array = data_log.get_column("x")
    
    # Plot 1: Temperature vs Distance
    plt.plot(x_array, data_log.get_column("T_L"))
    plt.title("Coolant Temperature Profile")
    plt.xlabel("Axial Distance (m)")
    plt.ylabel("Temperature (K)")
    
    # Plot 2: Mass Flow & Depletion Mechanics
    plt.plot(x_array, data_log.get_column("Gamma_L"), label="Total Flow")
    plt.title("Coolant Depletion Profile")
    
    # Plot 3: Evaporation vs Entrainment (Efficiency check)
    plt.plot(x_array, data_log.get_column("m_dot_Evap"), label="Evaporated (Useful)")
    plt.plot(x_array, data_log.get_column("m_dot_Ent"), label="Entrained (Wasted)")
    plt.legend()
    
    plt.show()