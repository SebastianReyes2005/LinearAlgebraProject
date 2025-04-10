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

spy_data = load_stock_data("SPY.csv")
bond_20_data = load_stock_data("TLT.csv")
# Assume the first command-line argument is the lower bound start date
start_limit = pd.Timestamp("2005")  # e.g., "2004-01-01"

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
        'Stock Gain Intraday': round(get_current_gain(bond_20_data, current_date) * 100, 5)
    })
gap_df = pd.DataFrame(gap_data)

def get_state(row, df):
    if gap_df.loc[row, 'Stock Gain Intraday'] >= 0:
        # up up
        if gap_df.loc[row, 'Stock Gap'] >= 0 and gap_df.loc[row, 'Bond Gap'] >= 0:  return 1
        # up down
        elif gap_df.loc[row, 'Stock Gap'] >= 0 and gap_df.loc[row, 'Bond Gap'] < 0: return 2
        # down up
        elif gap_df.loc[row, 'Stock Gap'] < 0 and gap_df.loc[row, 'Bond Gap'] >= 0: return 3
        # down down
        elif gap_df.loc[row, 'Stock Gap'] < 0 and gap_df.loc[row, 'Bond Gap'] < 0: return 4
    else:
        # up up
        if gap_df.loc[row, 'Stock Gap'] >= 0 and gap_df.loc[row, 'Bond Gap'] >= 0:  return 5
        # up down
        elif gap_df.loc[row, 'Stock Gap'] >= 0 and gap_df.loc[row, 'Bond Gap'] < 0: return 6
        # down up
        elif gap_df.loc[row, 'Stock Gap'] < 0 and gap_df.loc[row, 'Bond Gap'] >= 0: return 7
        # down down
        elif gap_df.loc[row, 'Stock Gap'] < 0 and gap_df.loc[row, 'Bond Gap'] < 0: return 8


a = np.zeros((8, 8), dtype=float)

for m in range(1, len(gap_data)):
    for column in range(1,9):
        current_count = 0
        if get_state(m-1, gap_df) == column:
            for row in range(1,9):
                if get_state(m, gap_df) == row:
                    a[row-1, column-1] += 1
print(a)
total_incidence = int(np.sum(a))
# print(total_incidence)
col_sums = np.sum(a, axis=0)
transition_matrix = np.zeros_like(a)
for col in range(8):
    if col_sums[col] > 0:  # avoid division by zero if no transitions from that state
        transition_matrix[:, col] = a[:, col] / col_sums[col]
print("\nNormalized Transition Matrix:")
print(transition_matrix)

# First, we will look at the gap probabilities in m=3
s_0 = [[0],[0],[0],[1],[0],[0],[0],[0]]
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
for i in range(1, 5):
    print(v[i-1]/(v[i-1]+v[(4+i)-1]))

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





