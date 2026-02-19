import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA

import torch
import torch.nn as nn
import torch.optim as optim

# ========= CONFIG =========
CSV_NAME    = "nn_3lock_yz_dataset.csv"  # dataset file name

EPOCHS      = 1000
BATCH_SIZE  = 8
PATIENCE    = 40        # early stopping patience
TRAIN_FRAC  = 0.85      # 85% train+val, 15% test
VAL_FRAC_OF_TRAIN = 0.1765  # inner val fraction (~15% of full data)
SPLIT_SEED  = 123       # controls data split only

LR          = 5e-4
WEIGHT_DECAY = 1e-4

N_PCS       = 8         # number of PCA components for Y (22 -> 8)
TOL_SWEEP   = [5.0, 10.0, 15.0, 20.0]  # NOW means percent thresholds
# ==========================


def find_csv_file(target_name: str) -> str:
    """
    Search the folder where this script is run (recursively)
    and return the full path to the CSV.
    """
    print("\n===== SEARCHING FOR DATASET =====")
    print(f"Looking for '{target_name}' starting at:\n{os.getcwd()}\n")

    for root, dirs, files in os.walk(os.getcwd()):
        if target_name in files:
            full_path = os.path.join(root, target_name)
            print(f"✔ Found CSV at: {full_path}\n")
            return full_path

    print("❌ Could not find the dataset file.")
    print("Here are some files I see in the current folder:")
    for f in os.listdir():
        print(" -", f)
    raise FileNotFoundError(f"\n'{target_name}' was not found anywhere under {os.getcwd()}\n")


def load_nn_dataset(csv_path: str):
    """
    Load nn_3lock_yz_dataset.csv and return (X_aug, Y_rel) arrays.

    X_aug: 15 inputs per sample
    Y_rel: 22 outputs (dy2,dz2,...,dy12,dz12)
    """
    df = pd.read_csv(csv_path)

    lock_cols = [f"L{i}" for i in range(1, 13)]
    y_cols    = [f"y{i}" for i in range(1, 13)]
    z_cols    = [f"z{i}" for i in range(1, 13)]

    for col in lock_cols + y_cols + z_cols:
        if col not in df.columns:
            raise ValueError(f"Required column '{col}' missing from dataset.")

    # Base lock bits
    X_lock = df[lock_cols].values.astype(np.float32)  # (N,12)
    N = X_lock.shape[0]

    # ---- Global features ----
    N_L = X_lock.sum(axis=1)

    indices = np.arange(1, 13, dtype=np.float32)
    eps = 1e-8
    mu = (X_lock * indices).sum(axis=1) / (N_L + eps)
    mu_norm = mu / 12.0

    dmax_list = []
    for row in X_lock:
        locked_idxs = np.where(row > 0.5)[0]
        if len(locked_idxs) == 0:
            dmax_list.append(0.0)
        else:
            dmax_list.append(float(locked_idxs.max() - locked_idxs.min()))
    dmax = np.array(dmax_list, dtype=np.float32)
    dmax_norm = dmax / 11.0

    global_feats = np.stack([N_L, mu_norm, dmax_norm], axis=1).astype(np.float32)

    # X_aug: shape (N,15)
    X_aug = np.concatenate([X_lock, global_feats], axis=1)

    # ---- Absolute Y ----
    Y_abs = df[y_cols + z_cols].values.astype(np.float32)  # (N,24)
    Y_abs_reshaped = Y_abs.reshape(N, 12, 2)

    ref = Y_abs_reshaped[:, 0:1, :]
    Y_rel_all = Y_abs_reshaped - ref
    Y_rel = Y_rel_all[:, 1:, :]  # drop cell 1
    Y_rel_flat = Y_rel.reshape(N, 22)

    return X_aug.astype(np.float32), Y_rel_flat.astype(np.float32)


class AuxeticNet(nn.Module):
    """2×8 MLP"""
    def __init__(self, in_dim=15, hidden_dims=[8, 8], out_dim=N_PCS):
        super().__init__()
        layers = []
        prev = in_dim
        for h in hidden_dims:
            layers.append(nn.Linear(prev, h))
            layers.append(nn.ReLU())
            prev = h
        layers.append(nn.Linear(prev, out_dim))
        self.layers = nn.Sequential(*layers)

    def forward(self, x):
        return self.layers(x)


# =====================================================================
#   PERCENTAGE ERROR VERSION (REPLACES OLD DISTANCE-BASED FUNCTIONS)
# =====================================================================

def compute_percentage_errors(Y_true, Y_pred):
    """
    %error = (||pred - true|| / ||true||) * 100
    If true vector is (0,0), denominator becomes 1e-8.
    """
    N = Y_true.shape[0]
    Yt = Y_true.reshape(N, 11, 2)
    Yp = Y_pred.reshape(N, 11, 2)

    true_norm = np.linalg.norm(Yt, axis=2) + 1e-8
    err_norm  = np.linalg.norm(Yp - Yt, axis=2)

    pct = (err_norm / true_norm) * 100.0
    return pct


def accuracy_at_percentage(pct_errors, max_pct):
    """
    Counts nodes with percentage error <= specified percent.
    """
    correct = (pct_errors <= max_pct).sum()
    total   = pct_errors.size
    return correct / total


# =====================================================================


def train_model(X_train, Z_train_n, X_val, Z_val_n, device):
    """
    Train one 2×8 model in PCA space.
    """
    Xtr = torch.tensor(X_train, device=device)
    Ztr = torch.tensor(Z_train_n, device=device)
    Xv  = torch.tensor(X_val, device=device)
    Zv  = torch.tensor(Z_val_n, device=device)

    model = AuxeticNet(in_dim=15, hidden_dims=[8, 8], out_dim=N_PCS).to(device)
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)

    best_val = float("inf")
    best_state = None
    no_imp = 0
    N_train = Xtr.shape[0]

    print("\n--- Training PCA model ---")

    for epoch in range(1, EPOCHS + 1):

        model.train()
        perm = torch.randperm(N_train)
        Xsh = Xtr[perm]
        Zsh = Ztr[perm]

        for i in range(0, N_train, BATCH_SIZE):
            xb = Xsh[i:i + BATCH_SIZE]
            zb = Zsh[i:i + BATCH_SIZE]

            optimizer.zero_grad()
            pred = model(xb)
            loss = criterion(pred, zb)
            loss.backward()
            optimizer.step()

        # Validation
        model.eval()
        with torch.no_grad():
            val_loss = criterion(model(Xv), Zv).item()

        if val_loss < best_val - 1e-6:
            best_val = val_loss
            best_state = model.state_dict()
            no_imp = 0
        else:
            no_imp += 1

        if epoch % 50 == 0 or epoch == 1:
            print(f"Epoch {epoch:4d} | Val Loss (PCA space): {val_loss:.4f}")

        if no_imp >= PATIENCE:
            print(f"Early stopping at epoch {epoch}, best val loss = {best_val:.4f}")
            break

    if best_state:
        model.load_state_dict(best_state)

    return model


def main():
    # ---------------- Load dataset ----------------
    path = find_csv_file(CSV_NAME)
    X_aug, Y_rel = load_nn_dataset(path)

    print("\nLoaded dataset with shapes:")
    print("X_aug:", X_aug.shape, " Y_rel:", Y_rel.shape)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print("\nUsing device:", device)

    # ---------------- Train/Test Split ----------------
    X_train_val, X_test, Y_train_val, Y_test = train_test_split(
        X_aug, Y_rel,
        test_size=1.0 - TRAIN_FRAC,
        random_state=SPLIT_SEED,
        shuffle=True
    )

    print("\nSplit sizes:")
    print("Train+Val:", X_train_val.shape[0])
    print("Test:     ", X_test.shape[0])

    # Inner split
    X_train, X_val, Y_train, Y_val = train_test_split(
        X_train_val, Y_train_val,
        test_size=VAL_FRAC_OF_TRAIN,
        random_state=SPLIT_SEED,
        shuffle=True
    )

    print("\nInner split:")
    print("Train:", X_train.shape[0])
    print("Val:  ", X_val.shape[0])

    # ---------------- PCA ----------------
    print(f"\nFitting PCA with {N_PCS} components...")
    pca = PCA(n_components=N_PCS)
    Z_train = pca.fit_transform(Y_train)
    Z_val   = pca.transform(Y_val)
    Z_test  = pca.transform(Y_test)

    print("Explained variance ratio:", pca.explained_variance_ratio_)
    print("Cumulative:", np.cumsum(pca.explained_variance_ratio_))

    # Normalize PCA coefficients
    Z_mean = Z_train.mean(axis=0, keepdims=True)
    Z_std  = Z_train.std(axis=0, keepdims=True) + 1e-8

    Ztr_n = (Z_train - Z_mean) / Z_std
    Zva_n = (Z_val   - Z_mean) / Z_std
    Zte_n = (Z_test  - Z_mean) / Z_std

    # ---------------- Train ----------------
    model = train_model(X_train, Ztr_n, X_val, Zva_n, device)

    # ---------------- Evaluate ----------------
    model.eval()
    Xte_t = torch.tensor(X_test, device=device)
    Zte_t = torch.tensor(Zte_n, device=device)

    with torch.no_grad():
        Zpred_n = model(Xte_t).cpu().numpy()

    # De-normalize
    Zpred = Zpred_n * Z_std + Z_mean

    # Back to 22-dim Y
    Y_pred = pca.inverse_transform(Zpred)
    Y_true = Y_test

    # Percentage errors
    pct_errors = compute_percentage_errors(Y_true, Y_pred)
    pct_flat = pct_errors.flatten()

    rmse_pct = np.sqrt(np.mean(pct_flat**2))
    mae_pct  = np.mean(np.abs(pct_flat))

    print("\n===== FINAL RESULTS (% ERROR METRIC) =====")
    print(f"RMSE (%): {rmse_pct:.2f}")
    print(f"MAE  (%): {mae_pct:.2f}")

    print("\nAccuracy vs percentage threshold:")
    for pct in TOL_SWEEP:
        acc = accuracy_at_percentage(pct_errors, pct)
        print(f"  <= {pct:4.1f}% error → {acc*100:6.2f}%")

    print()


if __name__ == "__main__":
    main()
