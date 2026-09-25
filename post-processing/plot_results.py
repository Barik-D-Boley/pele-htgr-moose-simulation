import pandas as pd
import matplotlib.pyplot as plt

def plot_temperature_and_stress():
  # Load MOOSE CSV output
  df = pd.read_csv('../moose/pele_core_sim_out.csv')
  
  fig, ax1 = plt.subplots(figsize=(8, 5))

  # Temperature plot
  color = 'tab:red'
  ax1.set_xlabel('Step / Time')
  ax1.set_ylabel('Max Temperature (K)', color=color)
  ax1.plot(df['time'], df['max_T'], color=color, linewidth=2, label='Max Temp')
  ax1.tick_params(axis='y', labelcolor=color)

  # Von Mises Stress plot on secondary axis
  ax2 = ax1.twinx()  
  color = 'tab:blue'
  ax2.set_ylabel('Max Von Mises Stress (Pa)', color=color)
  ax2.plot(df['time'], df['max_vonmises'], color=color, linewidth=2, linestyle='--', label='Max Stress')
  ax2.tick_params(axis='y', labelcolor=color)

  plt.title('HTGR RVE Thermal-Mechanical Response')
  fig.tight_layout()
  plt.savefig('rve_thermal_stress_summary.png', dpi=300)
  print("Plot saved to rve_thermal_stress_summary.png")

if __name__ == "__main__":
  plot_temperature_and_stress()