import click
import tqdm
from pathlib import Path
import tkinter as tk
from tkinter import filedialog

BAD_HEADERS = [b"ID", b"QUAL", b"FILTER", b"INFO", b"FORMAT"]


@click.group()
def cli():
    pass


def filter_meta_data_all_lines(all_lines):
    for i in range(len(all_lines)):
        if not all_lines[i].startswith(b"##"):
            break
    return all_lines[i:]


def filter_bad_snips(data_lines, filter_threshold):
    good_lines = []
    for line in data_lines:
        missing_percentage = (line.count(b"") / len(line)) * 100.0
        if missing_percentage <= filter_threshold:
            good_lines.append(line)
    return good_lines


def get_first_sample_index(header_line):
    return header_line.index(b"ALT") + 1
    for index, line in enumerate(header_line):
        try:
            int(line)
            return index
        except BaseException:
            continue
    raise ValueError("Failed to find sample index")


def filter_bad_samples(parsed_lines, filter_threshold):
    sample_index = get_first_sample_index(parsed_lines[0])
    indexes_to_remove = []
    for i in range(sample_index, len(parsed_lines[0])):
        missing_count = 0
        for l in parsed_lines[1:]:
            if l[i] == b"":
                missing_count += 1
        missing_percentage = (missing_count / len(parsed_lines)) * 100
        if missing_percentage > filter_threshold:
            indexes_to_remove.append(i)
    indexes_to_remove.sort(reverse=True)
    for l in parsed_lines:
        for i in indexes_to_remove:
            l.pop(i)


def filter_bad_headers(parsed_lines):
    indexes_to_remove = []
    header_line = parsed_lines[0]
    for bad in BAD_HEADERS:
        try:
            indexes_to_remove.append(header_line.index(bad))
        except ValueError:
            continue
    indexes_to_remove.sort(reverse=True)
    for l in parsed_lines:
        for i in indexes_to_remove:
            l.pop(i)


def replace_all(data_lines, pattern, replacement):
    for line in data_lines:
        for i, cell in enumerate(line):
            if cell.startswith(pattern):
                line[i] = replacement


def add_counters(data_lines, values_to_count):
    for line in data_lines:
        for v in values_to_count:
            line.append(str(line.count(v)).encode())


def add_multiplication(data_lines):
    for line in data_lines:
        a, b, c = line[-3], line[-2], line[-1]
        a = int(a) + 1
        b = int(b) + 1
        c = int(c) + 1
        line.append(str(a * b * c).encode())


def add_minor(data_lines):
    for line in data_lines:
        a, b, c = float(line[-4]), float(line[-3]), float(line[-2])
        total = (c + (b / 2.0)) / float(a + b + c)
        line.append(f"{total:.2f}".encode())


def group_linkage(data_lines, linkage_limit=350):
    pos_index = data_lines[0].index(b"POS")
    first_sample_index = get_first_sample_index(data_lines[0])
    groups = [data_lines[0]]
    cur = []
    for i, line in enumerate(data_lines[1:]):
        if len(cur) == 0:
            cur.append(line)
        elif abs(int(line[pos_index]) - int(cur[0][pos_index])) < linkage_limit:
            cur.append(line)
        else:
            groups.append(flatten_linkage(first_sample_index, cur))
            cur = [line]
    if len(cur) > 0:
        groups.append(flatten_linkage(first_sample_index, cur))

    data_lines.clear()
    data_lines.extend(groups)


def flatten_linkage(first_sample_index, linkage_lines):
    flat = linkage_lines[0][:first_sample_index]
    for i in range(first_sample_index, len(linkage_lines[0])):
        counters = [0, 0, 0]
        total_count = 0
        for line in linkage_lines:
            if len(line[i]) > 0:
                counters[int(line[i]) - 1] += 1
        total = sum(counters)
        for i, val in filter(lambda x: x[1] > total / 2.0, enumerate(counters)):
            flat.append(str(i + 1).encode())
            break
        else:
            flat.append(b"")

    return flat


FILTER_THRESHOLD = 6

REPLACEMENTS = (
    (b"./.", b""),
    (b"0/0:", b"1"),
    (b"0/1:", b"2"),
    (b"1/1:", b"3"),
    (b'"0/0:', b"1"),
    (b'"0/1:', b"2"),
    (b'"1/1:', b"3"),
    (b"NC_057849.1", b"1"),
)


def t_diff(lines, src, dst, filt, sample_index):
    name = lines[0][src]
    res = []
    for i in range(1, len(lines)):
        if (
            len(lines[i][src]) == 0
            or len(lines[i][dst]) == 0
            or lines[i][src] == b"--"
            or lines[i][dst] == b"--"
            or lines[i][dst] == b"./."
            or lines[i][src] == b"./."
        ):
            continue
        res.append(abs(int(lines[i][src]) - int(lines[i][dst])))
    if len(res) < filt:
        return ""
    return "{:.5f}".format(sum(res) / float(len(res)))


PEARSONS_DISTANCE_MAP = {
    (1, 1): 1,
    (2, 2): 1,
    (3, 3): 1,
    (1, 2): 0.7071067,
    (2, 1): 0.7071067,
    (3, 2): 0.7071067,
    (2, 3): 0.7071067,
    (1, 3): 0,
    (3, 1): 0,
}


def pearsons_diff(lines, src, dst, filt, sample_index):
    name = lines[0][src]
    res = []
    for i in range(1, len(lines)):
        if (
            len(lines[i][src]) == 0
            or len(lines[i][dst]) == 0
            or lines[i][src] == b"--"
            or lines[i][dst] == b"--"
            or lines[i][dst] == b"./."
            or lines[i][src] == b"./."
        ):
            continue
        res.append(PEARSONS_DISTANCE_MAP[(int(lines[i][src]), int(lines[i][dst]))])
    if len(res) < filt:
        return ""
    return "{:.5f}".format(sum(res) / float(len(res)))


MORISITA_DISTANCE_MAP = {
    (1, 1): 1,
    (3, 3): 1,
    (2, 2): 1,
    (1, 2): 2/3,
    (2, 1): 2/3,
    (3, 2): 2/3,
    (2, 3): 2/3,
    (1, 3): 0,
    (3, 1): 0,
}


def morisitas_diff(lines, src, dst, filt, sample_index):
    name = lines[0][src]
    res = []
    for i in range(1, len(lines)):
        if (
            len(lines[i][src]) == 0
            or len(lines[i][dst]) == 0
            or lines[i][src] == b"--"
            or lines[i][dst] == b"--"
            or lines[i][dst] == b"./."
            or lines[i][src] == b"./."
        ):
            continue
        res.append(MORISITA_DISTANCE_MAP[(int(lines[i][src]), int(lines[i][dst]))])
    if len(res) < filt:
        return ""
    return "{:.5f}".format(sum(res) / float(len(res)))

def t_shared(lines, src, dst, filt, sample_index):
    name = lines[0][src]
    total_snips = len(lines) - 1
    count = 0
    for i in range(1, len(lines)):
        if (
            len(lines[i][src]) == 0
            or len(lines[i][dst]) == 0
            or lines[i][src] == b"--"
            or lines[i][dst] == b"--"
            or lines[i][dst] == b"./."
            or lines[i][src] == b"./."
        ):
            continue
        count += 1
    return "{:.2f}".format(count / total_snips)


def calculate_distances(lines, diff_function=t_diff):
    res = {}
    sample_index = get_first_sample_index(lines[0])
    for src in tqdm.tqdm(range(sample_index, len(lines[0]))):
        res[lines[0][src]] = []
        for dst in range(sample_index, len(lines[0])):
            if src != dst:
                res[lines[0][src]].append(diff_function(lines, src, dst, 0, sample_index))
            else:
                res[lines[0][src]].append("")
    return res


def calculate_shared_snps(lines):
    res = {}
    sample_index = get_first_sample_index(lines[0])
    for src in tqdm.tqdm(range(sample_index, len(lines[0]))):
        res[lines[0][src]] = []
        for dst in range(sample_index, len(lines[0])):
            if src != dst:
                res[lines[0][src]].append(t_shared(lines, src, dst, 0, sample_index))
            else:
                res[lines[0][src]].append("")
    return res


@cli.command()
@click.argument("snip_filter", type=float)
@click.argument("sample_filter", type=float)
@click.option("--input_vcf", "-i", type=click.File("rb"), default=None)
@click.option("--limit", "-l", type=int, default=None)
def filter_lines(snip_filter, sample_filter, input_vcf, limit):
    """
    Filter and process vcf file

    input_vcf : The vcf file to process.

    snip_filter: The allowed percentage of missing data before filtering snip

    sample_filter: The allowed percentage of missing data before filtering sample
    """
    if input_vcf is None:
        root = tk.Tk()
        root.withdraw()
        file_path = filedialog.askopenfilename()
        input_vcf = click.File("rb")(file_path)
    all_data = input_vcf.read()
    all_lines = all_data.split(b"\n")
    all_lines = filter_meta_data_all_lines(all_lines)
    parsed_data = [l.split(b"\t") for l in all_lines[:limit]]
    # Remove empty last line
    parsed_data.pop()
    filter_bad_headers(parsed_data)
    print(f"Total lines: {len(parsed_data)}")
    for pat, val in REPLACEMENTS:
        replace_all(parsed_data, pat, val)

    if click.confirm("linkage?"):
        group_linkage(parsed_data)
        print(f"After linkage grouping: {len(parsed_data)}")
    else:
        print("No linkage!")

    good_lines = filter_bad_snips(parsed_data, snip_filter)
    print(f"After Snip filter: {len(good_lines)}")
    filter_bad_samples(good_lines, sample_filter)

    header_line = good_lines.pop(0)
    add_counters(good_lines, (b"1", b"2", b"3"))
    add_multiplication(good_lines)
    add_minor(good_lines)
    header_line.extend([b"111", b"222", b"333", b"123", b"mnr"])
    good_lines.insert(0, header_line)
    good_lines = [b"\t".join(l) for l in good_lines]

    out_path = Path(input_vcf.name)
    out_path = out_path.with_name(
        out_path.name.replace(".vcf", "_parsed.vcf").replace(".txt", "_parsed.txt")
    )
    out_path.write_bytes(b"\n".join(good_lines))


@cli.command()
@click.option("--input_vcf", "-i", type=click.File("rb"), default=None)
@click.option("--limit", "-l", type=int, default=None)
@click.option("--pearson", "-p", is_flag=True, default=False)
@click.option("--morisita", "-m", is_flag=True, default=False)
def distances(input_vcf, limit, pearson, morisita):
    """
    Calculates distance table for each sample by sample

    input_vcf : The vcf file to process.
    """
    if input_vcf is None:
        root = tk.Tk()
        root.withdraw()
        file_path = filedialog.askopenfilename()
        input_vcf = click.File("rb")(file_path)
    all_data = input_vcf.read()
    all_lines = all_data.split(b"\n")
    parsed_data = [l.split(b"\t")[:-5] for l in all_lines[:limit]]
    parsed_data.pop()
    if pearson and morisita:
        print("Can't do both Pearson's and Morisita's at the same time, please only choose one (-p/-m)")
    elif morisita:
        print("Doing Morisita's™ distance")
        distances = calculate_distances(parsed_data, morisitas_diff)
    elif pearson:
        print("Doing Pearson's™ distance")
        distances = calculate_distances(parsed_data, pearsons_diff)
    else:
        print("Doing Yaron's™ distance")
        distances = calculate_distances(parsed_data)
    out_path = Path(input_vcf.name)
    if morisita:
        out_path = out_path.with_name(out_path.name + "_morisita_distances.csv")
    elif pearson:
        out_path = out_path.with_name(out_path.name + "_pearson_distances.csv")
    else:
        out_path = out_path.with_name(out_path.name + "_distances.csv")
    out = out_path.open("wb")
    first_line = ",".join(name.decode() for name in distances.keys())
    first_line = "," + first_line
    first_line += "\r\n"
    out.write(first_line.encode())

    for name, v in distances.items():
        line = name.decode() + ","
        line += ",".join(str(x) for x in v)
        line += "\r\n"
        out.write(line.encode())
    out.close()


@cli.command()
@click.option("--input_vcf", "-i", type=click.File("rb"), default=None)
@click.option("--limit", "-l", type=int, default=None)
def shared_snps(input_vcf, limit):
    """
    Calculates shared snips table for each sample by sample

    input_vcf : The vcf file to process.
    """
    if input_vcf is None:
        root = tk.Tk()
        root.withdraw()
        file_path = filedialog.askopenfilename()
        input_vcf = click.File("rb")(file_path)
    all_data = input_vcf.read()
    all_lines = all_data.split(b"\n")
    parsed_data = [l.split(b"\t")[:-5] for l in all_lines[:limit]]
    parsed_data.pop()
    shared_snps = calculate_shared_snps(parsed_data)
    out_path = Path(input_vcf.name)
    out_path = out_path.with_name(out_path.name + "_shared_snps.csv")
    out = out_path.open("wb")
    first_line = ",".join(name.decode() for name in shared_snps.keys())
    first_line = "," + first_line
    first_line += "\r\n"
    out.write(first_line.encode())

    for name, v in shared_snps.items():
        line = name.decode() + ","
        line += ",".join(str(x) for x in v)
        line += "\r\n"
        out.write(line.encode())
    out.close()


if __name__ == "__main__":
    cli()
