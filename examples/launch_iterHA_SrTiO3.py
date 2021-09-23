import click
from phonopy import Phonopy
from phonopy.interface.vasp import read_vasp_from_strings
from aiida.manage.configuration import load_profile
from aiida.orm import QueryBuilder, Int, Float, Bool, Str, load_node
from aiida.plugins import DataFactory, WorkflowFactory
from aiida.engine import run, submit
from aiida_phonopy.common.utils import (
    phonopy_atoms_to_structure,
    phonopy_atoms_from_structure,
)
from aiida_phonopy.workflows.iter_ha import _extract_dataset_from_db, _create_dataset

load_profile()


def launch_aiida():

    Dict = DataFactory("dict")
    unitcell_str = """ Sr Ti O
   1.0
     3.9050000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000    3.9050000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000    3.9050000000000000
 Sr Ti O
   1   1   3
Direct
   0.0000000000000000  0.0000000000000000  0.0000000000000000
   0.5000000000000000  0.5000000000000000  0.5000000000000000
   0.5000000000000000  0.0000000000000000  0.5000000000000000
   0.5000000000000000  0.5000000000000000  0.0000000000000000
   0.0000000000000000  0.5000000000000000  0.5000000000000000"""

    cell = read_vasp_from_strings(unitcell_str)
    structure = phonopy_atoms_to_structure(cell)

    base_incar_dict = {
        "PREC": "Accurate",
        "IBRION": -1,
        "EDIFF": 1e-8,
        "NELMIN": 5,
        "NELM": 100,
        "ENCUT": 500,
        "IALGO": 38,
        "ISMEAR": 0,
        "SIGMA": 0.01,
        "GGA": "PS",
        "LREAL": False,
        "lcharg": False,
        "lwave": False,
    }

    base_config = {
        "code_string": "vasp544mpi@nancy",
        "potential_family": "PBE.54",
        "potential_mapping": {"O": "O", "Ti": "Ti_pv", "Sr": "Sr_sv"},
        "options": {
            "resources": {
                "num_machines": 1,
                "parallel_env": "mpi*",
                "tot_num_mpiprocs": 24,
            },
            "max_wallclock_seconds": 3600 * 10,
        },
    }
    base_parser_settings = {
        "add_energies": True,
        "add_forces": True,
        "add_stress": True,
    }
    forces_config = base_config.copy()
    forces_config.update(
        {
            "kpoints_mesh": [4, 4, 4],  # k-point density,
            "parser_settings": base_parser_settings,
            "parameters": base_incar_dict,
        }
    )
    nac_config = base_config.copy()
    nac_parser_settings = {"add_born_charges": True, "add_dielectrics": True}
    nac_parser_settings.update(base_parser_settings)
    nac_incar_dict = {"lepsilon": True}
    nac_incar_dict.update(base_incar_dict)
    nac_config.update(
        {
            "kpoints_mesh": [8, 8, 8],  # k-point density,
            "parser_settings": nac_parser_settings,
            "parameters": nac_incar_dict,
        }
    )

    # PhononPhonopy = WorkflowFactory('phonopy.phonopy')
    # builder = PhononPhonopy.get_builder()
    PhonopyIterHA = WorkflowFactory("phonopy.iter_ha")
    builder = PhonopyIterHA.get_builder()
    builder.structure = structure
    builder.calculator_settings = Dict(
        dict={"forces": forces_config, "nac": nac_config}
    )
    builder.run_phonopy = Bool(False)
    builder.remote_phonopy = Bool(True)
    builder.code_string = Str("phonopy@nancy")
    builder.phonon_settings = Dict(
        dict={
            "mesh": 50.0,
            "supercell_matrix": [2, 2, 2],
            "distance": 0.01,
            "is_nac": True,
            "fc_calculator": "alm",
        }
    )
    builder.symmetry_tolerance = Float(1e-5)
    builder.options = Dict(dict=base_config["options"])
    builder.metadata.label = "SrTiO3 iterative phonon 2x2x2 1000K test"
    builder.metadata.description = "SrTiO3 iterative phonon 2x2x2 1000K test"
    builder.max_iteration = Int(50)
    builder.number_of_snapshots = Int(40)
    builder.temperature = Float(1000.0)
    builder.number_of_steps_for_fitting = Int(20)
    builder.include_ratio = Float(0.99)
    # builder.initial_nodes = Dict(
    #     dict={'nodes':  [86164, 86936, 87708, 88480, 89252,
    #                      90024, 90796, 91568, 92340, 93112]})

    # Chose how to run the calculation
    run_by_deamon = True
    if not run_by_deamon:
        result = run(builder)
        print(result)
    else:
        future = submit(builder)
        print(future)
        print("Running workchain with pk={}".format(future.pk))


def get_nac_params(pk_nac):
    n_nac = load_node(pk_nac)
    if "nac_params" in n_nac.outputs:
        borns = n_nac.outputs.nac_params.get_array("born_charges")
        epsilon = n_nac.outputs.nac_params.get_array("epsilon")
        nac_params = {"born": borns, "factor": 14.399652, "dielectric": epsilon}
        return nac_params
    else:
        return None


def get_phonon(pk, pk_nac):
    n = load_node(pk)
    unitcell = phonopy_atoms_from_structure(n.inputs.structure)
    smat = n.outputs.phonon_setting_info["supercell_matrix"]
    ph = Phonopy(unitcell, supercell_matrix=smat, primitive_matrix="auto")
    force_sets = n.outputs.force_sets.get_array("force_sets")
    dataset = n.outputs.phonon_setting_info["displacement_dataset"]
    ph.dataset = dataset
    ph.forces = force_sets
    ph.nac_params = get_nac_params(pk_nac)

    return ph


def find_latest_uuid():
    IterHarmonicApprox = WorkflowFactory("phonopy.iter_ha")
    qb = QueryBuilder()
    qb.append(IterHarmonicApprox)
    qb.order_by({IterHarmonicApprox: {"ctime": "desc"}})
    qb.first()
    return qb.first()[0].uuid


def search_pk(uuid):
    """uuid can be pk."""
    IterHarmonicApprox = WorkflowFactory("phonopy.iter_ha")
    qb = QueryBuilder()
    qb.append(IterHarmonicApprox, tag="iter_ha", filters={"uuid": {"==": uuid}})
    PhonopyWorkChain = WorkflowFactory("phonopy.phonopy")
    qb.append(PhonopyWorkChain, with_incoming="iter_ha")
    qb.order_by({PhonopyWorkChain: {"ctime": "asc"}})
    pks = [n[0].pk for n in qb.all() if n[0].is_finished_ok]

    return pks


def get_num_prev(uuid):
    n = load_node(uuid)
    return n.inputs.number_of_steps_for_fitting.value


def get_initial_nodes(uuid):
    n = load_node(uuid)
    if "initial_nodes" in n.inputs:
        return n.inputs["initial_nodes"]["nodes"]
    else:
        return None


def get_include_ratio(uuid):
    n = load_node(uuid)
    if "include_ratio" in n.inputs:
        return n.inputs.include_ratio.value
    else:
        return None


def get_temperature(uuid):
    n = load_node(uuid)
    return n.inputs.temperature.value


def bunch_phonons(pks, pk_nac, max_items=None, include_ratio=None, linear_decay=True):
    nodes = [load_node(pk) for pk in pks]
    displacements, forces, energies = _extract_dataset_from_db(
        [n.outputs.force_sets for n in nodes],
        [n.outputs.phonon_setting_info for n in nodes],
    )
    d, f, e, idx = _create_dataset(
        displacements,
        forces,
        energies,
        max_items,
        include_ratio,
        linear_decay=linear_decay,
    )
    unitcell = phonopy_atoms_from_structure(nodes[0].inputs.structure)
    smat = nodes[0].outputs.phonon_setting_info["supercell_matrix"]
    ph = Phonopy(unitcell, supercell_matrix=smat, primitive_matrix="auto")
    ph.dataset = {"displacements": d, "forces": f}
    ph.nac_params = get_nac_params(pk_nac)
    return ph


def bunch_band_phonopy(
    pks, pk_nac, max_items=None, include_ratio=None, linear_decay=True
):
    ph = bunch_phonons(
        pks,
        pk_nac,
        max_items=max_items,
        include_ratio=include_ratio,
        linear_decay=linear_decay,
    )
    print(
        "max_items: %d, include_ratio: %f, linear_decay: %s, snapshots: %d"
        % (max_items, include_ratio, linear_decay, ph.dataset["displacements"].shape[0])
    )
    ph.produce_force_constants(fc_calculator="alm")
    pks_str = get_pks_str(pks)
    filename = "band-%s.yaml" % pks_str
    ph.auto_band_structure(write_yaml=True, filename=filename)
    # filename = "band-%s.hdf5" % pks_str
    # ph.auto_band_structure()
    # ph.write_hdf5_band_structure(filename=filename)
    print("%s was made for PKs=%s." % (filename, pks_str))


def band_phonopy(pk, pk_nac):
    ph = get_phonon(pk, pk_nac)
    ph.produce_force_constants(fc_calculator="alm")
    filename = "band-%d.yaml" % pk
    ph.auto_band_structure(write_yaml=True, filename=filename)
    print("%s was made for PK=%d." % (filename, pk))


def get_bunch_tp_phonopy(
    pks, pk_nac, temperature, max_items=None, include_ratio=None, linear_decay=True
):
    ph = bunch_phonons(
        pks,
        pk_nac,
        max_items=max_items,
        include_ratio=include_ratio,
        linear_decay=linear_decay,
    )
    ph.produce_force_constants(fc_calculator="alm")
    ph.run_mesh(mesh=100.0, shift=[0.5, 0.5, 0.5])
    ph.run_thermal_properties(
        temperatures=[
            temperature,
        ]
    )
    free_energy = ph.get_thermal_properties_dict()["free_energy"][0]
    return free_energy


def get_tp_phonopy(pk, pk_nac, temperature):
    ph = get_phonon(pk, pk_nac)
    ph.produce_force_constants(fc_calculator="alm")
    ph.run_mesh(mesh=100.0, shift=[0.5, 0.5, 0.5])
    ph.run_thermal_properties(
        temperatures=[
            temperature,
        ]
    )
    free_energy = ph.get_thermal_properties_dict()["free_energy"][0]
    return free_energy


def bunch_dump_phonopy(
    pks, pk_nac, max_items=None, include_ratio=None, linear_decay=True
):
    ph = bunch_phonons(
        pks,
        pk_nac,
        max_items=max_items,
        include_ratio=include_ratio,
        linear_decay=linear_decay,
    )
    settings = {
        "force_sets": True,
        "displacements": True,
        "force_constants": False,
        "born_effective_charge": True,
        "dielectric_constant": True,
    }
    pks_str = get_pks_str(pks)
    filename = "phonopy_params-%s.yaml" % pks_str
    ph.save(filename=filename, settings=settings)
    print("%s was made for PKs=%s." % (filename, pks_str))


def dump_phonopy(pk, pk_nac):
    ph = get_phonon(pk, pk_nac)
    settings = {
        "force_sets": True,
        "displacements": True,
        "force_constants": False,
        "born_effective_charge": True,
        "dielectric_constant": True,
    }
    filename = "phonopy_params-%d.yaml" % pk
    ph.save(filename=filename, settings=settings)
    print("%s was made for PK=%d." % (filename, pk))


def get_pks_str(pks):
    if len(pks) > 4:
        pks_str = "%d-%dpks-%d" % (pks[0], len(pks) - 2, pks[-1])
    else:
        pks_str = "-".join(["%d" % pk for pk in pks])
    return pks_str


def get_pks_list(pks, num_prev, initial_pks=None, num_tail=None):
    if initial_pks is None:
        pks_list = [
            [
                pks[0],
            ],
        ]
        _pks = pks[1:]
    else:
        pks_list = []
        _pks = initial_pks + pks

    for i in range(len(_pks)):
        pks_appended = _pks[: (i + 1)]
        if len(pks_appended) > num_prev:
            pks_appended = pks_appended[-num_prev:]
        if initial_pks is None:
            pks_list.append(pks_appended)
        else:
            if i > len(initial_pks) - 2:
                pks_list.append(pks_appended)

    if num_tail is None:
        return pks_list
    else:
        if len(pks_list) <= num_tail:
            return pks_list
        else:
            return pks_list[-num_tail:]


@click.command()
@click.argument("pks", nargs=-1)
@click.option("--launch", default=False, is_flag=True)
@click.option("--pk", default=None, type=int)
@click.option("--uuid", default=None)
@click.option("--dump", "bunch_dump", default=False, is_flag=True)
@click.option("--band", "bunch_band", default=False, is_flag=True)
@click.option("--fe", "bunch_free_energy", default=False, is_flag=True)
@click.option("--list-pks", "list_pks", default=False, is_flag=True)
@click.option("--show-uuid", "show_uuid", default=False, is_flag=True)
@click.option("--pk-nac", "pk_nac", default=None, type=int)
@click.option("--num-prev", "num_prev", default=None, type=int)
@click.option("--tail", "num_tail", default=None, type=int)
@click.option("--include-ratio", "include_ratio", default=None, type=float)
@click.option("--linear-decay", "linear_decay", default=False, is_flag=True)
def main(
    launch,
    pk,
    uuid,
    pks,
    bunch_dump,
    bunch_band,
    bunch_free_energy,
    list_pks,
    show_uuid,
    pk_nac,
    num_prev,
    num_tail,
    include_ratio,
    linear_decay,
):
    """

    Launching the calculation, PK of IterHA is shown. This is good to
    rembember, e.g., keeping it in AiiDA Groups. Then we can quickly view how
    many iterations have already done by

    % verdi process show PK_of_IterHA

    ...
    Called        PK  Type
    --------  ------  ----------------
    CALL      117655  WorkChainNode
    CALL      117742  CalcFunctionNode
    CALL      117747  WorkChainNode
    CALL      118514  CalcFunctionNode
    CALL      118519  WorkChainNode
    CALL      119286  CalcFunctionNode
    ...

    where each WorkChainNode node corresponds to one iteration.
    Unspecified builder.initial_nodes, the first one is the 0K phonon
    calculation. With is_nac=True, this should only contain nac_params.
    This PK is used to specify --pk-nac in this script.

    Usage
    -----
    The initial calculation can be launched by --launch option.

    % python launch_phonon_SrTiO3.py --launch

    The following options are used to analyize the (on-going) results.
    --pk or --uuid is used to specify which IterHA calculation will be
    analyzed. When omitting it, the last IterHA calculation's uuid is
    searched, but it is always better to specify it since this search
    can fail.

    To list PKs of PhonopyWorkChain calculations

    % python launch_phonon_SrTiO3.py --pk [PK] --list-pks

    """

    if launch:
        launch_aiida()
    else:
        if pk is not None:
            _uuid = load_node(pk).uuid
            print("PK: %d" % pk)
        elif uuid is None:
            _uuid = find_latest_uuid()
        else:
            _uuid = uuid

        print("UUID: %s" % _uuid)

        if pk_nac is None:
            _pk_nac = search_pk(_uuid)[0]
        else:
            _pk_nac = pk_nac

        if num_prev is None:
            _num_prev = get_num_prev(_uuid)
        else:
            _num_prev = num_prev

        if include_ratio is None:
            _include_ratio = get_include_ratio(_uuid)
        else:
            _include_ratio = include_ratio

        if pks and not list_pks:
            _pks = [int(pk) for pk in pks]
        else:
            _pks = search_pk(_uuid)

        initial_nodes = get_initial_nodes(_uuid)

        if show_uuid:  # Show latest IterHarmonicApprox node uuid
            print(find_latest_uuid())
        elif list_pks:  # List PKs of PhonopyWorkChain
            print(
                "All PhonopyWorkChain nodes: ",
                ",".join(["%d" % pk for pk in search_pk(_uuid)]),
            )
            if len(_pks) > 0:
                print(
                    "PhonopyWorkChain nodes used to generate random "
                    "displacements at each iteration:"
                )
                pks_list = get_pks_list(
                    _pks,
                    get_num_prev(_uuid),
                    initial_pks=initial_nodes,
                    num_tail=num_tail,
                )
                if initial_nodes is not None:
                    print("Initial nodes: %s" % initial_nodes)
                else:
                    print("%d (0K)" % pks_list.pop(0)[0])
                for pk_set in pks_list:
                    print(",".join(["%d" % pk for pk in pk_set]))
        elif bunch_dump and _pks:
            for pk_set in get_pks_list(
                _pks, _num_prev, initial_pks=initial_nodes, num_tail=num_tail
            ):
                bunch_dump_phonopy(
                    pk_set,
                    _pk_nac,
                    max_items=_num_prev,
                    include_ratio=_include_ratio,
                    linear_decay=linear_decay,
                )
        elif bunch_band and _pks:
            for pk_set in get_pks_list(
                _pks, _num_prev, initial_pks=initial_nodes, num_tail=num_tail
            ):
                bunch_band_phonopy(
                    pk_set,
                    _pk_nac,
                    max_items=_num_prev,
                    include_ratio=_include_ratio,
                    linear_decay=linear_decay,
                )
        elif bunch_free_energy and _pks:
            temperature = get_temperature(_uuid)
            for i, pk_set in enumerate(
                get_pks_list(
                    _pks, _num_prev, initial_pks=initial_nodes, num_tail=num_tail
                )
            ):
                fe = get_bunch_tp_phonopy(
                    pk_set,
                    _pk_nac,
                    temperature,
                    max_items=_num_prev,
                    include_ratio=_include_ratio,
                    linear_decay=linear_decay,
                )
                print("%d %f" % (i + 1, fe))


if __name__ == "__main__":
    main()
