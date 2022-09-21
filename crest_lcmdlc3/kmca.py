import subprocess
import sys
from glob import glob
from io import StringIO

import conda.cli
import matplotlib

matplotlib.use("Agg")
import importlib

import matplotlib.pyplot as plt
import numpy as np
import scipy.cluster.hierarchy
from scipy.spatial.distance import euclidean, squareform
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


def conda_install(package):
    proc = subprocess.run(
        ["conda", "install"] + [package],
        text=True,
        capture_output=True,
    )
    return proc.stdout


rmsd_loader = importlib.util.find_spec("rmsd")
if rmsd_loader is not None:
    from rmsd import main as calculate_rmsd
else:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "rmsd"])
    from rmsd import main as calculate_rmsd
ase_loader = importlib.util.find_spec("ase")
if ase_loader is not None:
    import ase
    from ase.io import read, write, iread
else:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "ase"])
    import ase
    from ase.io import read, write, iread
tblite_loader = importlib.util.find_spec("tblite")
if tblite_loader is not None:
    import tblite
    from tblite.ase import TBLite
    calculator = TBLite
else :
    from ase.calculators.lj import LennardJones
    calculator = LennardJones()


# MAIN SETUP OPTIONS #
fileroot = sys.argv[1].split(".")[0]
n_clusters = 10
dirname = f"{fileroot}_selected"
method = "mix"
commonpath = os.getcwd()
save = True
######################

filenames = []
for i,mol in enumerate(iread(sys.argv[1])):
    mol.write(f"{fileroot}_conformer_{i}.xyz")
    filenames.append(f"{fileroot}_conformer_{i}.xyz")


print(
    "WARNING! This script seriously assumes that all structures are equally sorted in terms of atom indices. If this is not the case, contact ruben.laplazasolanas@epfl.ch"
)
# print(f"Filenames are: {filenames} with base directory {commonpath}.")


class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio  # free up some memory
        sys.stdout = self._stdout


def plot_dendrogram(m: np.ndarray, label: str):
    # Start clustering based off rmsds
    print(
        f"Max pairwise {label}: {np.max(m)} in whatever units the argument was passed."
    )
    assert np.all(m - m.T < 1e-6)
    reduced_distances = squareform(m)
    linkage = scipy.cluster.hierarchy.linkage(reduced_distances, method="average")
    plt.title(f"{label} Average linkage hierarchical clustering")
    dn = scipy.cluster.hierarchy.dendrogram(
        linkage, no_labels=True, count_sort="descendent"
    )
    plt.savefig(f"{label}_dendrogram.png")


def do_kmeans(n_clusters: int, m: np.ndarray, filenames: list):
    km = KMeans(n_clusters=n_clusters, n_init=50)
    scaler = StandardScaler()
    m = scaler.fit_transform(m)
    cm = km.fit_predict(m)
    u, c = np.unique(cm, return_counts=True)
    closest_pt_idx = []
    print(f"Unique clusters found: {u} \nWith counts: {c} \nAdding up to {np.sum(c)}")
    for iclust in range(km.n_clusters):
        # get all points assigned to each cluster:
        cluster_pts = m[km.labels_ == iclust]
        # get all indices of points assigned to this cluster:
        cluster_pts_indices = np.where(km.labels_ == iclust)[0]

        cluster_cen = km.cluster_centers_[iclust]
        min_idx = np.argmin(
            [euclidean(m[idx], cluster_cen) for idx in cluster_pts_indices]
        )

        # Testing:
        print(
            f"Closest index of point to cluster {iclust} center: {filenames[cluster_pts_indices[min_idx]]}"
        )
        closest_pt_idx.append(cluster_pts_indices[min_idx])
    return [filenames[cpi] for cpi in closest_pt_idx]


def compute_rmsd_m(filenames: list):
    distances = np.zeros((len(filenames), len(filenames)))
    for i, _ in enumerate(filenames):
        for j in range(i, len(filenames)):
            if i == j:
                distances[i, j] = 0
            else:
                with Capturing() as output:
                    calculate_rmsd([filenames[i], filenames[j], "-nh"])
                distances[i, j] = output[0]
                distances[j, i] = output[0]
    return distances


def compute_rele_m(filenames: list):
    distances = np.zeros((len(filenames), len(filenames)))
    es = np.zeros((len(filenames)))
    for i, _ in enumerate(filenames):
        mol = read(filenames[i])
        # coords = mol.get_positions()
        mol.calc = calculator
        es[i] = mol.get_potential_energy() * 23.0609
    for i, _ in enumerate(filenames):
        for j in range(i, len(filenames)):
            if i == j:
                distances[i, j] = 0
            else:
                distances[i, j] = np.abs(es[i] - es[j])
                distances[j, i] = np.abs(es[i] - es[j])
    return distances


if __name__ == "__main__":
    if method == "rmsd":
        print("Working with rmsds in angstrom.")
        m = compute_rmsd_m(filenames)
        plot_dendrogram(m, label="rmsd")
        selected = do_kmeans(n_clusters, m, filenames)
        if save:
            destination = f"{dirname}_rmsd"
            print(f"Saving selected conformers to {destination}.")
            if not os.path.exists(f"{destination}"):
                os.mkdir(f"{destination}")
            for filename in selected:
                os.popen(f"cp {filename} {destination}")

    if method == "rele":
        print("Working with relative energies in kcal/mol.")
        m = compute_rele_m(filenames)
        plot_dendrogram(m, label="rele")
        selected = do_kmeans(n_clusters, m, filenames)
        if save:
            destination = f"{dirname}_rele"
            print(f"Saving selected conformers to {destination}.")
            if not os.path.exists(f"{destination}"):
                os.mkdir(f"{destination}")
            for filename in selected:
                os.popen(f"cp {filename} {destination}")

    if method == "mix":
        print("Working with hybrid similarity metric.")
        m1 = compute_rmsd_m(filenames)
        m2 = compute_rele_m(filenames)
        m1 = np.abs(m1) / np.max(m1)
        m2 = np.abs(m2) / np.max(m2)
        m = (np.sqrt(m1) + np.sqrt(m2)) / 2
        selected = do_kmeans(n_clusters, m, filenames)
        if save:
            destination = f"{dirname}_mix/"
            print(f"Saving selected conformers to {destination}.")
            if not os.path.exists(f"{destination}"):
                os.mkdir(f"{destination}")
            for filename in selected:
                os.popen(f"cp {filename} {destination}")
