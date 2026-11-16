import pytest

from scbatchflow.simulate import SimConfig, simulate


@pytest.fixture(scope="session")
def sim_dir(tmp_path_factory):
    out = tmp_path_factory.mktemp("sim")
    simulate(out, SimConfig(cells_per_sample=40, empty_droplets=15, mean_umis=120,
                            n_background_genes=40, seed=7))
    return out
