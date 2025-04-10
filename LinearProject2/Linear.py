import numpy as np
import pandas as pd
import datetime as dt
import plotly.graph_objects as go


def load_stock_data(file_path):
    """
    Load CSV with columns: Date, High, Low, Adj Close.
    """
    df = pd.read_csv(file_path, parse_dates=['Date'])
    df.sort_values('Date', inplace=True)
    df.set_index('Date', inplace=True)
    return df

def get_current_price(df, current_date):
    """Return the latest Adj Close on or before current_date."""
    data = df.loc[df.index <= current_date]
    if data.empty:
        return np.nan
    return data.iloc[-1]['Close']

def get_current_gain(df, current_date):
    data = df.loc[df.index <= current_date]
    if data.empty:
        return np.nan
    return (data.iloc[-1]['Close'] - data.iloc[-1]['Open']) / data.iloc[-1]['Open']

def get_gap(df, current_date):
    data = df.loc[df.index <= current_date]
    if len(data) < 2 or data.empty:
        return np.nan
    return (data.iloc[-1]['Open'] - data.iloc[-2]['Close'])/ data.iloc[-2]['Close']

import os
base_dir = os.path.dirname(os.path.abspath(__file__))
spy_path = os.path.join(base_dir, "SPY.csv")
tlt_path = os.path.join(base_dir, "TLT.csv")

spy_data = load_stock_data(spy_path)
bond_20_data = load_stock_data(tlt_path)
# Assume the first command-line argument is the lower bound start date
start_limit = pd.Timestamp("2004-12-01")  # e.g., "2004-01-01"


# Now, for each DataFrame, keep only the data on or after start_limit.
# bond_20_data = bond_20_data.loc[bond_20_data.index >= start_limit]
spy_data = spy_data.loc[spy_data.index >= start_limit]
gap_data = []

for i in range(len(spy_data)):
    current_date = spy_data.index[i]
    gap_data.append({
        'Date': current_date,
        'Stock Gap': round(get_gap(spy_data, current_date) * 100, 5),
        'Bond Gap': round(get_gap(bond_20_data, current_date) * 100, 5),
        'Stock Gain Intraday': round(get_current_gain(spy_data, current_date) * 100, 5)
    })
gap_df = pd.DataFrame(gap_data)


# def get_state(row, df):
#     # up up
#     if gap_df.loc[row, 'Stock Gap'] >= 0 and gap_df.loc[row, 'Bond Gap'] >= 0:  return 1
#     # up down
#     elif gap_df.loc[row, 'Stock Gap'] >= 0 and gap_df.loc[row, 'Bond Gap'] < 0: return 2
#     # down up
#     elif gap_df.loc[row, 'Stock Gap'] < 0 and gap_df.loc[row, 'Bond Gap'] >= 0: return 3
#     # down down
#     elif gap_df.loc[row, 'Stock Gap'] < 0 and gap_df.loc[row, 'Bond Gap'] < 0: return 4

def get_state(row, df):
    z_stock = (gap_df.loc[row, 'Stock Gap'] - np.mean(gap_df['Stock Gap']))/np.std(gap_df['Stock Gap'])
    z_bond = (gap_df.loc[row, 'Bond Gap'] - np.mean(gap_df['Bond Gap']))/np.std(gap_df['Bond Gap'])
    # up up
    if z_stock >= 0 and z_bond >= 0:  return 1
    # up down
    elif z_stock >= 0 and z_bond < 0: return 2
    # down up
    elif z_stock < 0 and z_bond >= 0: return 3
    # down down
    elif z_stock < 0 and z_bond < 0: return 4





    

a = np.zeros((4, 4), dtype=float)

for m in range(1, len(gap_data)):
    for column in range(1,5):
        current_count = 0
        if get_state(m-1, gap_df) == column:
            for row in range(1,5):
                if get_state(m, gap_df) == row:
                    a[row-1, column-1] += 1
print (a)
total_incidence = int(np.sum(a))
# print(total_incidence)
col_sums = np.sum(a, axis=0)
transition_matrix = np.zeros_like(a)
for col in range(4):
    if col_sums[col] > 0:  # avoid division by zero if no transitions from that state
        transition_matrix[:, col] = a[:, col] / col_sums[col]
print("\nNormalized Transition Matrix:")
print(transition_matrix)

# First, we will look at the gap probabilities in m=3
s_0 = [[1],[0],[0],[0]]
s_1 = np.matmul(transition_matrix,s_0)
print(s_1)
s_2 = np.matmul(transition_matrix,s_1)
print(s_2)




# Algebraic Method of Steady Vector
I = np.eye(transition_matrix.shape[0])
# Construct the homogeneous system: (I - p)v = 0
A = I - transition_matrix
# To pick a unique solution, we replace one of the equations with the constraint:
# x1 + x2 + x3 + x4 = 1
# One approach is to augment the system by appending the constraint row.
A_aug = np.vstack([A, np.ones(transition_matrix.shape[0])])
b_aug = np.zeros(transition_matrix.shape[0] + 1)
b_aug[-1] = 1  # the right-hand side of the sum constraint is 1
# Solve the augmented system using least squares:
v, residuals, rank, s = np.linalg.lstsq(A_aug, b_aug, rcond=None)
print("Steady-state vector v:")
print(v)

# Large Sample Method of Steady State (by raising the transition matrix to a very high power)
# Using np.linalg.matrix_power ensures numerical stability.
long_run = np.linalg.matrix_power(transition_matrix, 10000)
print("Long-run Transition Matrix:")
print(long_run)


            

# trade_df = gap_df[(gap_df['Stock Gap'] >= 0) &
#                   ((gap_df['TLT Gap']  >= 0) | (gap_df['TLT Gap'].isna()))]


# # trade_df = gap_df[(gap_df['Stock Gap'] <= -3.0)]

# # print(gap_df.to_string(index=True))
# print(trade_df.to_string(index=True))

def round_matrix_sf(matrix, sf):
  """Rounds all elements in a matrix to a specified number of significant figures.

  Args:
    matrix: A NumPy array representing the matrix.
    sf: The number of significant figures to round to.

  Returns:
    A new NumPy array with the rounded values.
  """
  rounded_matrix = np.array([[round_sf(element, sf) for element in row] for row in matrix])
  return rounded_matrix

def round_sf(x, sf):
  """Rounds a single number to a specified number of significant figures."""
  if x == 0:
    return 0
  magnitude = np.floor(np.log10(abs(x)))
  factor = 10**(sf - 1 - magnitude)
  return round(x * factor) / factor

# Step 1: Extract 'Stock Gap' values into a raw array
raw = np.zeros((40, 1), dtype=float)
for i in range(40):
    raw[i, 0] = gap_df.loc[i + 1, 'Stock Gap']

# Step 2: Compute moving average (days 1–20, 2–21, ..., 21–40)
x = np.zeros((10, 1), dtype=float)
for i in range(10):
    x[i, 0] = np.average(raw[i:i + 20, 0])

# Step 3: Set polynomial degree
degree = 9  # should be <= len(x) - 1 for best fit

# Step 4: Construct matrix M for polynomial fitting
m = np.zeros((len(x), degree), dtype=float)
for column in range(degree):
    for row in range(len(x)):
        m[row, column] = (row + 1) ** column
print("\nWe get M")
print(round_matrix_sf(m,5))

# Step 5: Solve for polynomial coefficients Z using least squares
print("\nWe get M^T")
m_T = m.T
print(round_matrix_sf(m_T,5))
print("\nWe get (M^T)M")
m_T_m = np.matmul(m_T, m)
print(round_matrix_sf(m_T_m,5))
print("\nWe get Z=(((M^T)M)^-1)X")
z = np.matmul(np.matmul(np.linalg.inv(m_T_m), m_T), x)
print(round_matrix_sf(z,5))

# Step 6: Reconstruct fitted values (in-sample check)
print("\n--- Reconstructed Values vs Original x ---")
result = []
for day in range(len(x)):
    fitted = sum(z[i, 0] * (day + 1) ** i for i in range(degree))
    result.append({'Day': day + 1, 'Original': round(x[day, 0], 6), 'Fitted': round(fitted, 6)})
result_df = pd.DataFrame(result)
print(result_df.to_string(index=False))