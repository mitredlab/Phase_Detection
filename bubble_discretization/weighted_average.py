import pandas as pd

# Given table data
data = {
    'Frequency': [184, 110, 59, 31, 11, 7, 3, 2, 2, 1],
    'PRE Area (%)': [-0.50635, 0.18401, 0.02811, 0.01039, -0.00787, 0.00266, 0.00266, 0.00643, 0.00643, 0.00643],
    'PRE Perimeter (%)': [59.79746, 23.98589, 18.06246, 15.62654, 14.33143, 12.92963, 12.92963, 11.99070, 11.99070, 11.99070],
    'ME Area (µm^2)': [-1.59e-12, 5.60e-12, 2.41e-12, 1.76e-12, -2.21e-12, 1.56e-12, 1.56e-12, 8.08e-12, 8.08e-12, 8.08e-12],
    'ME Perimeter (µm)': [3.76e-05, 4.69e-05, 5.93e-05, 7.20e-05, 8.50e-05, 1.11e-04, 1.11e-04, 1.51e-04, 1.51e-04, 1.51e-04]
}

# Creating a DataFrame
df = pd.DataFrame(data)

# Calculating weighted averages
weighted_avg_PRE_Area = (df['Frequency'] * df['PRE Area (%)']).sum() / df['Frequency'].sum()
weighted_avg_PRE_Perimeter = (df['Frequency'] * df['PRE Perimeter (%)']).sum() / df['Frequency'].sum()
weighted_avg_ME_Area = (df['Frequency'] * df['ME Area (µm^2)']).sum() / df['Frequency'].sum()
weighted_avg_ME_Perimeter = (df['Frequency'] * df['ME Perimeter (µm)']).sum() / df['Frequency'].sum()

weighted_avg_PRE_Area, weighted_avg_PRE_Perimeter, weighted_avg_ME_Area, weighted_avg_ME_Perimeter
