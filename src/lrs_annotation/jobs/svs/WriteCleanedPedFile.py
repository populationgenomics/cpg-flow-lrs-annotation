import re
from cpg_utils import Path, to_path
from cpg_utils.config import ConfigError, reference_path

PED_FAMILY_ID_REGEX = re.compile(r'(^[A-Za-z0-9_]+$)')

def get_references(keys: list[str | dict[str, str]]) -> dict[str, str | list[str]]:
    """
    Dict of WDL inputs with reference file paths.
    """
    res: dict[str, str | list[str]] = {}
    for key in keys:
        # Keys can be maps (e.g. {'MakeCohortVcf.cytobands': 'cytoband'})
        if isinstance(key, dict):
            key, ref_d_key = next(iter(key.items()))  # noqa: PLW2901
        else:
            ref_d_key = key
        # e.g. GATKSVPipelineBatch.rmsk -> rmsk
        ref_d_key = ref_d_key.split('.')[-1]
        try:
            res[key] = reference_path(f'gatk_sv/{ref_d_key}')
        except KeyError:
            res[key] = reference_path(f'broad/{ref_d_key}')
        except ConfigError:
            res[key] = reference_path(f'broad/{ref_d_key}')

    return res

def clean_ped_family_ids(ped_line: str) -> str:
    """
    Takes each line in the pedigree and cleans it up
    If the family ID already conforms to expectations, no action
    If the family ID fails, replace all non-alphanumeric/non-underscore
    characters with underscores

    >>> clean_ped_family_ids('family1\tchild1\t0\t0\t1\t0\\n')
    'family1\tchild1\t0\t0\t1\t0\\n'
    >>> clean_ped_family_ids('family-1-dirty\tchild1\t0\t0\t1\t0\\n')
    'family_1_dirty\tchild1\t0\t0\t1\t0\\n'

    Args:
        ped_line (str): line from the pedigree file, unsplit

    Returns:
        the same line with a transformed family id
    """

    split_line = ped_line.rstrip().split('\t')

    if re.match(PED_FAMILY_ID_REGEX, split_line[0]):
        return ped_line

    # if the family id is not valid, replace failing characters with underscores
    split_line[0] = re.sub(r'[^A-Za-z0-9_]', '_', split_line[0])

    # return the rebuilt string, with a newline at the end
    return '\t'.join(split_line) + '\n'


def make_clean_combined_ped(ped_file_path: Path, output_path: Path) -> Path:
    """
    Create cohort + ref panel PED.
    Concatenating all samples across all datasets with ref panel

    There are restrictions on valid characters in PED file, so we clean up the family IDs
    to conform to the regex: ^[A-Za-z0-9_]+$.
    """
    conf_ped_path = get_references(['ped_file'])['ped_file']
    assert isinstance(conf_ped_path, str)
    with output_path.open('w') as out:
        with ped_file_path.open() as f:
            # layer of family ID cleaning
            for line in f:
                out.write(clean_ped_family_ids(line))
        # The ref panel PED doesn't have any header, so can safely concatenate:
        with to_path(conf_ped_path).open() as f:
            out.write(f.read())
    return output_path
