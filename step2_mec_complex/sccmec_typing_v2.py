import os
import argparse
from argparse import RawTextHelpFormatter
from executor import ExternalCommand
import pandas as pd
import json
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Align
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from datetime import datetime


def path_convert(p):
    p = p.replace("\\", "/")
    return p


def execute(cmd, directory=os.getcwd(), capture=False, capture_stderr=False,
            stdout_file=None, stderr_file=None, silent=True):
    """A simple wrapper around executor."""
    # from executor import ExternalCommand
    command = ExternalCommand(
        cmd, directory=directory, capture=capture,
        capture_stderr=capture_stderr, stdout_file=stdout_file,
        stderr_file=stderr_file, silent=silent
    )
    command.start()


def exist_mk(p):
    makedir(p)


def makedir(p):
    if os.path.exists(p):
        # print("文件已存在")
        pass
    else:
        execute(f"mkdir {p}")


def make_blast_db_genome(input_genome, path_db_dir):
    genome_name = input_genome.split(".")[0]
    title = f"blastdb_{genome_name}"
    path_db = path_convert(os.path.join(path_db_dir, title))
    path_title = path_convert(os.path.join(path_db_dir, title))
    path_genome = path_convert(os.path.join(args.path_input_genome, input_genome))
    execute(f"makeblastdb -in {path_genome} -dbtype nucl -title {path_title} -parse_seqids -out {path_db}")
    return path_title

def make_blast_db_scc(input_scc, path_db_dir):
    scc_name = input_scc.split("/")[-1].split(".")[0]
    title = f"blastdb_{scc_name}"
    path_db = path_convert(os.path.join(path_db_dir, title))
    path_title = path_convert(os.path.join(path_db_dir, title))
    execute(f"makeblastdb -in {input_scc} -dbtype nucl -title {path_title} -parse_seqids -out {path_db}")
    return path_title

def find_alleles(path_queries, path_db, path_blast_result, genome, threshold):
    genome_name = genome.split(".")[0]
    df_queries = pd.DataFrame()
    for query in os.listdir(path_queries):
        query_name = query.replace(".fa", "")
        path_query = path_convert(os.path.join(path_queries, query))
        if "scc" in path_db:
            file_output = blast_query_scc(genome_name, query_name, path_query, path_db, path_blast_result)
        else:
            file_output = blast_query(genome_name, query_name, path_query, path_db, path_blast_result)
        result_line = convert_blast_result(genome_name, query_name, file_output)

        df_query = pd.DataFrame(result_line)
        if df_queries.shape[0] == 0:
            df_queries = df_query
        else:
            df_queries = pd.concat([df_queries, df_query], axis=0)
    df_queries["coverage"] = df_queries["coverage"].astype(float)
    # coverage的阈值筛选
    condition = (df_queries['contig'] != "0") & (df_queries['coverage'] > threshold)
    df_filtered = df_queries[condition]
    df_result_final = pd.DataFrame()
    for contig, group in df_filtered.groupby('contig'):
        df_result_contig = pd.DataFrame()
        group = group.reset_index(drop=True)  # 确保分组内索引唯一
        # 比较所有行（使用 iloc 访问）
        for index1, row_i in group.iterrows():
            df_row = pd.DataFrame(row_i).T
            if df_result_contig.shape[0] == 0:
                df_result_contig = df_row
            else:
                judge_new = True
                for index2, row_j in df_result_contig.iterrows():
                    # 提取数值型参数
                    a_start = min(row_i['hit_from'], row_i['hit_to'])
                    a_end = max(row_i['hit_from'], row_i['hit_to'])
                    a_contig = row_i["contig"]
                    b_start = min(row_j['hit_from'], row_j['hit_to'])
                    b_end = max(row_j['hit_from'], row_j['hit_to'])
                    b_contig = row_j["contig"]
                    if is_same_region(a_start, a_end, b_start, b_end, 200, a_contig, b_contig):
                        judge_new = False
                        if row_i["coverage"] > row_j["coverage"]:
                            # 替换
                            df_result_contig.loc[index2] = row_i
                        else:
                            pass
                if judge_new:
                    df_result_contig = pd.concat([df_result_contig, df_row], axis=0)
        if df_result_contig.shape[0] == 0:
            df_result_contig = pd.DataFrame({"genome": [genome_name], "query_name": [f"no_{query_name}"], "identity":[0], "coverage":[0],
                                             "hit_from":[0], "hit_to":[0], "query_from": [0], "query_to": [0], "strand":[0], "contig":[0]})
        if df_result_final.shape[0] == 0:
            df_result_final = df_result_contig
        else:
            df_result_final = pd.concat([df_result_final, df_result_contig], axis=0)
    if df_result_final.shape[0] == 0:
        df_result_final = pd.DataFrame(
            {"genome": [genome_name], "query_name": [f"no_{query_name}"], "identity": [0], "coverage": [0],
             "hit_from": [0], "hit_to": [0], "query_from": [0], "query_to": [0],"strand": [0], "contig": [0]})
    df_result_final = df_result_final.reset_index(drop=True)
    return df_result_final


def is_same_region(a_start, a_end, b_start, b_end, threshold, a_contig, b_contig):
    if a_contig == b_contig:
        # 完全包含
        if ((a_start <= b_start) and (a_end >= b_end)) or ((b_start <= a_start) and (b_end >= a_end)):
            return True
        # 相交且非重叠部分总长度 < 100
        overlap_start = max(a_start, b_start)
        overlap_end = min(a_end, b_end)
        if overlap_start <= overlap_end:  # 存在重叠
            total_length = (a_end - a_start) + (b_end - b_start)
            non_overlap = total_length - 2 * (overlap_end - overlap_start)
            if non_overlap < threshold:
                return True
            else:
                return False
        return False
    else:
        return False


def blast_query(genome_name, query_name, path_query, path_db, path_blast_result):
    file_output = path_convert(os.path.join(path_blast_result, f"{genome_name}_{query_name}.json"))
    execute(f"blastn -query {path_query} -out {file_output} -outfmt 15 -db {path_db} -evalue 0.001")
    # execute(f"blastn -query {path_query} -task blastn-short -out {file_output} -outfmt 15 -db {path_db} -evalue 0.001 -word_size 7 -reward 1 -penalty -1 -gapopen 2 -gapextend 1")

    return file_output

def blast_query_scc(genome_name, query_name, path_query, path_db, path_blast_result):
    file_output = path_convert(os.path.join(path_blast_result, f"{genome_name}_{query_name}.json"))
    # execute(f"blastn -query {path_query} -task blastn-short -out {file_output} -outfmt 15 -db {path_db} -evalue 0.001 -word_size 7 -reward 1 -penalty -1 -gapopen 2 -gapextend 1 -soft_masking false")
    execute(f"blastn -query {path_query} -out {file_output} -outfmt 15 -db {path_db} -evalue 0.001")
    return file_output
def convert_blast_result(genome_name, query_name, file_json):
    file = open(file_json)
    f = json.load(file)
    file.close()

    result_output = dict()
    result_output["genome"] = []
    result_output["query_name"] = []
    result_output["identity"] = []
    result_output["coverage"] = []
    result_output["hit_from"] = []
    result_output["hit_to"] = []
    result_output["query_from"] = []
    result_output["query_to"] = []
    result_output["strand"] = []
    result_output["contig"] = []

    querys = f["BlastOutput2"]
    for query in querys:
        results = query["report"]["results"]["search"]
        query_len = results["query_len"]
        hits = results["hits"]
        if hits:
            hsp_num = 0
            for hit in hits:
                hit_id = hit["description"][0]["id"]
                for hsp in hit["hsps"]:
                    result_output["genome"].append(genome_name)
                    result_output["query_name"].append(query_name)
                    result_output["identity"].append("{:.2f}".format(hsp["identity"] / hsp["align_len"]))
                    result_output["coverage"].append("{:.2f}".format(hsp["identity"] / query_len))
                    result_output["hit_from"].append(hsp["hit_from"])
                    result_output["hit_to"].append(hsp["hit_to"])

                    result_output["query_from"].append(hsp["query_from"])
                    result_output["query_to"].append(hsp["query_to"])

                    result_output["strand"].append(hsp["hit_strand"])
                    result_output["contig"].append(hit_id)
        else:
            result_output["genome"].append(genome_name)
            result_output["query_name"].append(query_name)
            result_output["identity"].append("{:.2f}".format(0.0))
            result_output["coverage"].append("{:.2f}".format(0.0))
            result_output["hit_from"].append("0")
            result_output["hit_to"].append("0")
            result_output["query_from"].append("0")
            result_output["query_to"].append("0")
            result_output["strand"].append("0")
            result_output["contig"].append("0")
    return result_output


def extract_from_genome(input_genome, query_name, df, start_shift, end_shift, path_scc_save):
    genome_name = input_genome.split(".")[0]
    genome_path = path_convert(os.path.join(args.path_input_genome, input_genome))
    records = list(SeqIO.parse(genome_path, "fasta"))
    scc_path = []
    for i, row in df.iterrows():
        strand = row["strand"]
        if strand == "Plus":
            start = row["hit_from"] - start_shift
            end = row["hit_to"] + end_shift
        else:
            start = row["hit_to"] - end_shift
            end = row["hit_from"] + start_shift
        if start < 0:
            start = 0
        records_saved = []
        for record in records:
            if record.id == row["contig"]:
                length = len(record.seq)
                if end > length:
                    end = length
                if strand == "Plus":
                    seq_new = record.seq[start: end]
                    record_new = SeqRecord(seq_new, id=f"{genome_name}_{query_name}_{i + 1}", description="")
                else:
                    seq_new = record.seq[start: end].reverse_complement()
                    record_new = SeqRecord(seq_new, id=f"{genome_name}_{query_name}_{i+1}_reversed", description="")
                records_saved.append(record_new)

        file_saved_path = path_convert(os.path.join(path_scc_save, f"{genome_name}_{query_name}_{i+1}.fa"))

        with open(file_saved_path, "w") as f:
            SeqIO.write(records_saved, f, "fasta")
        scc_path.append(file_saved_path)

    return scc_path


def find_previous_row(df, col, word):
    df = df.reset_index(drop=True)
    previous_row = pd.Series(
            {"genome": "0", "query_name": word, "identity": 0, "coverage": 0,
             "hit_from": 0, "hit_to": 0, "strand": 0, "contig": 0})
    # 找到query_name包含"mec"的行
    mec_row = df[df[col].str.contains(word)]
    # 获取mec行的索引
    if not mec_row.empty:
        mec_row_index = mec_row.index[0]
        # 检查是否是第一行
        if mec_row_index > 0:
            # 获取上一行
            previous_row = df.iloc[mec_row_index - 1]
    return previous_row



def judge_mec_complex(df_mec_scc, df_is_scc, df_meci_scc, df_meci_genome):
    id = df_mec_scc.iloc[0]["contig"]
    mec_name = df_mec_scc.iloc[0]["query_name"]
    type_mec_complex = "Unknown"
    note_complex = " "
    if mec_name == "mecC":
        type_mec_complex = "E"
    else:
        df_meci = pd.concat([df_meci_scc, df_meci_genome], axis=0)
        condition = (df_meci["query_name"] == "mecI") & (df_meci["coverage"] >= 0.95)
        df_meci_filted = df_meci[condition]
        df_mec_is = pd.concat([df_mec_scc, df_is_scc], axis=0).sort_values(by=["hit_from"],
                                                                           ascending=[True]).reset_index(drop=True)
        condition = (df_mec_is["hit_from"] < 4010) & (df_mec_is["coverage"] > 0)
        df_mec_is_filted = df_mec_is[condition].copy()
        previous_row = find_previous_row(df_mec_is_filted, "query_name", "mec")
        row_name, row_strand, row_converage = previous_row["query_name"], previous_row["strand"], previous_row["coverage"]
        # 有完整mecI
        if df_meci_filted.shape[0] >= 1:
            if row_strand != 0 and row_name == "IS431":
                if row_strand == "Minus":
                    type_mec_complex = "A3"
                else:
                    type_mec_complex = "A2"
                note_complex = f"{row_name}_{row_strand}_{row_converage}"
            else:
                type_mec_complex = "A1"
                note_complex = "no_IS"
        # 没有完整mecI，其他型
        else:
            if row_name == "IS1272":
                type_mec_complex = "B"
                note_complex = f"{row_name}_{row_strand}_{row_converage}"
            else:
                if row_name == "IS431":
                    if row_strand == "Plus":
                        type_mec_complex = "C1"
                    else:
                        type_mec_complex = "C2"
                    note_complex = f"{row_name}_{row_strand}_{row_converage}"
                else:

                    note_complex = previous_row["hit_from"]
                    if previous_row["hit_from"] < 10:
                        type_mec_complex = "Untypable"
                    else:
                        type_mec_complex = "Untypable or D"

    return id, type_mec_complex, note_complex

def df_result_emerge(mec_id, type_mec_complex, note_complex, df_mec_alleles_scc, df_meci_alleles_scc, df_meci_alleles_genome, df_is_alleles_scc, df_ccr_alleles_scc, df_ccr_alleles_genome):
    result = {"genome": [], "scc_id": [], "mec_complex": [], "note_complex":[], "mec_coverage":[],
              "mecI_scc":[], "mecI_genome":[], "mecR_scc":[], "mecR_genome":[],
              "IS_left":[], "IS_right":[], "ccr_scc":[], "ccr_genome":[]}
    genome_id = mec_id.split("_scc")[0]
    result["genome"].append(genome_id)
    result["scc_id"].append(mec_id)
    result["mec_complex"].append(type_mec_complex)
    result["note_complex"].append(note_complex)
    result["mec_coverage"].append(extract_gene_coverage(df_mec_alleles_scc, "mec"))
    result["mecI_scc"].append(extract_gene_coverage(df_meci_alleles_scc, "mecI"))
    result["mecI_genome"].append(extract_gene_coverage(df_meci_alleles_genome, "mecI"))
    result["mecR_scc"].append(extract_gene_coverage(df_meci_alleles_scc, "mecR"))
    result["mecR_genome"].append(extract_gene_coverage(df_meci_alleles_genome, "mecR"))
    if isinstance(note_complex, str):
        if note_complex.startswith("IS"):
            result["IS_left"].append(note_complex)
        else:
            result["IS_left"].append("no")
    else:
        result["IS_left"].append("no")
    result["IS_right"].append(is431_right(df_mec_alleles_scc, df_is_alleles_scc))
    result["ccr_scc"].append(extract_gene_coverage(df_ccr_alleles_scc, "ccr"))
    result["ccr_genome"].append(extract_gene_coverage(df_ccr_alleles_genome, "ccr"))
    df_res = pd.DataFrame(result)
    return df_res


def is431_right(df1, df2):
    df = pd.concat([df1, df2], axis=0).sort_values(by=["hit_from"],ascending=[True]).reset_index(drop=True)
    next_row = "empty"
    # 找到query_name包含"mec"的行
    mec_row = df[df["query_name"].str.contains("mec")]
    # 获取mec行的索引
    if not mec_row.empty:
        mec_row_index = mec_row.index[0]
        # 检查是否是最后一行
        if mec_row_index < (df.shape[0] -1):
            # 获取下一行
            next_row = df.iloc[mec_row_index + 1]
    if isinstance(next_row, str):
        content = "empty"
    else:
        gene_name = next_row["query_name"]
        row_strand = next_row["strand"]
        row_coverage = next_row["coverage"]
        content = f"{gene_name}_{row_strand}_{row_coverage}"
    return content


def extract_gene_coverage(df, gene_name):
    df = df[df["coverage"] > 0]
    df_filtered = df[df["query_name"].str.contains(gene_name)]
    content = ""
    if df_filtered.shape[0] > 0:
        gene_name_ser = df_filtered["query_name"]
        gene_coverage_ser = df_filtered["coverage"]
        for index, value in gene_name_ser.items():
            cover = gene_coverage_ser[index]
            content += f"{value}_{cover} "
    else:
        content = "No"
    return content

# def extract_coverage(df, gene_name):
#     print(gene_name)
#     print("*"* 40)
#     print(df)
#     df_filtered = df[df["query_name"].str.contains(gene_name)]
#     coverage = ""
#     print("*"* 40)
#     print(df_filtered)
#     if df_filtered.shape[0] == 0:
#         coverage = "0"
#     else:
#         coverage_series = df_filtered["query_name"]
#         for i, v in coverage_series.items():
#             coverage += v
#     print(coverage)
#     return coverage

def complete_seq_is(seq_new, row, orientation):
    query_name = row["query_name"]
    file_input = path_convert(os.path.join(args.scheme, f"IS/{query_name}.fa"))
    record = SeqIO.read(file_input, "fasta")

    if orientation == "left":
        if row["strand"] == "Plus":
            seq_queshi = record.seq[0: row["query_from"]]
        else:
            seq_queshi = record.seq[row["query_to"]:].reverse_complement()
        seq_new = seq_queshi + seq_new
    else:
        if row["strand"] == "Plus":
            seq_queshi = record.seq[row["query_to"]:]
        else:
            seq_queshi = record.seq[0: row["query_from"]:].reverse_complement()
        seq_new = seq_new + seq_queshi
    return seq_new

def extract_from_scc(file_input_seq, df_mec, df_ugpq, df_is, path_mec_complex_saved, path_hvr_save, path_except_hvr_save):
    record = SeqIO.read(file_input_seq, "fasta")
    records_mec_complex_saved = []
    records_HVR_saved = []
    records_e_HVR_saved = []
    df_all = pd.concat([df_mec, df_ugpq, df_is], axis=0).sort_values(by=["hit_from"],ascending=[True])
    df_all[["identity", "coverage"]] = df_all[["identity", "coverage"]].astype("float")
    df_all[["hit_from", "hit_to"]] = df_all[["hit_from", "hit_to"]].astype("int")
    condition = (df_all["coverage"] > 0) & (df_all["identity"] > 0.9)
    df_all = df_all[condition].reset_index(drop=True)

    # 找到ugpQ的下一个IS
    next_row = "empty"
    ugp_row = df_all[df_all["query_name"].str.contains("ugpQ")]
    # 获取mec行的索引
    if not ugp_row.empty:
        mec_row_index = ugp_row.index[0]
        # 检查是否是最后一行
        if mec_row_index < (df_all.shape[0] - 1):
            # 获取下一行
            next_row = df_all.iloc[mec_row_index + 1]
            if next_row["coverage"] < 0.6:
                if (mec_row_index + 2) <= (df_all.shape[0] - 1):
                    next_row = df_all.iloc[mec_row_index + 2]

    # 找到mecA的上一个IS
    previous_row = find_previous_row(df_all, "query_name", "mec")
    # 抽取mec complex
    complete = "complete"
    if isinstance(next_row, str):
        complex_end = len(record.seq)
        complete = "incomplete"
    else:
        complex_end = max(next_row["hit_from"], next_row["hit_to"])

    if previous_row["contig"] != 0:
        complex_start = min(previous_row["hit_from"], previous_row["hit_to"])
    else:
        complex_start = 0
        complete = "incomplete"

    seq_new = record.seq[complex_start: complex_end]
    if complete == "complete":
        seq_new = complete_seq_is(seq_new, previous_row, "left")
        seq_new = complete_seq_is(seq_new, next_row, "right")

    mec_id = file_input_seq.split("/")[-1].replace(".fa", "")
    new_file_name = f"{complete}_{mec_id}_mec_complex"
    record_new = SeqRecord(seq_new, id=new_file_name, description="")
    records_mec_complex_saved.append(record_new)
    file_saved_path = path_convert(os.path.join(path_mec_complex_saved, f"{new_file_name}.fa"))
    with open(file_saved_path, "w") as f:
        SeqIO.write(records_mec_complex_saved, f, "fasta")
    # hvr区域
    complete = "complete"
    if isinstance(next_row, str):
        hvr_end = len(record.seq)
        complete = "incomplete"
    else:
        hvr_end = max(next_row["hit_from"], next_row["hit_to"])

    if not ugp_row.empty:
        ugp_row = ugp_row.iloc[0]
        hvr_start = min(ugp_row["hit_from"], ugp_row["hit_to"])
    else:
        hvr_start = 0

    seq_new = record.seq[hvr_start: hvr_end]
    if complete == "complete":
        seq_new = complete_seq_is(seq_new, next_row, "right")

    mec_id = file_input_seq.split("/")[-1].replace(".fa", "")
    new_file_name = f"{complete}_{mec_id}_HVR"
    record_new = SeqRecord(seq_new, id=new_file_name, description="")
    records_HVR_saved.append(record_new)
    file_saved_path = path_convert(os.path.join(path_hvr_save, f"{new_file_name}.fa"))
    with open(file_saved_path, "w") as f:
        SeqIO.write(records_HVR_saved, f, "fasta")

    # e_hvr区域
    complete = "complete"
    if previous_row["contig"] != 0:
        e_hvr_start = min(previous_row["hit_from"], previous_row["hit_to"])
    else:
        e_hvr_start = 0
        complete = "incomplete"

    if not ugp_row.empty:
        e_hvr_end = min(ugp_row["hit_from"], ugp_row["hit_to"])
    else:
        e_hvr_end = len(record.seq)
        complete = "incomplete"

    seq_new = record.seq[e_hvr_start: e_hvr_end]
    if complete == "complete":
        seq_new = complete_seq_is(seq_new, previous_row, "left")

    mec_id = file_input_seq.split("/")[-1].replace(".fa", "")
    new_file_name = f"{complete}_{mec_id}_e_hvr"
    record_new = SeqRecord(seq_new, id=new_file_name, description="")
    records_e_HVR_saved.append(record_new)
    file_saved_path = path_convert(os.path.join(path_except_hvr_save, f"{new_file_name}.fa"))
    with open(file_saved_path, "w") as f:
        SeqIO.write(records_e_HVR_saved, f, "fasta")


# def judge_ccr_complex(df_ccr_scc, df_ccr_genome):
#     print(df_ccr_scc)
#     print(df_ccr_genome)
#     pass



def main():
    # 建立临时文件夹和输出文件夹
    # 结果文件夹
    path_result_dir = path_convert(os.path.join(args.output_dir, "result"))
    exist_mk(path_result_dir)
    # blast_db保存
    path_db_dir = path_convert(os.path.join(path_result_dir, "1.db"))
    exist_mk(path_db_dir)

    # blast结果保存
    path_blast_result = path_convert(os.path.join(path_result_dir, "2.blast_result"))
    exist_mk(path_blast_result)
    # scc seq保存
    path_scc_save = path_convert(os.path.join(path_result_dir, "3.scc_seq"))
    exist_mk(path_scc_save)
    # mecA复合体 seq保存
    path_mec_complex_save = path_convert(os.path.join(path_result_dir, "4.mec_complex_seq"))
    exist_mk(path_mec_complex_save)
    # HVR seq保存
    path_hvr_save = path_convert(os.path.join(path_result_dir, "5.HVR_seq"))
    exist_mk(path_hvr_save)
    # 保守 seq保存
    path_except_hvr_save = path_convert(os.path.join(path_result_dir, "6.except_HVR_seq"))
    exist_mk(path_except_hvr_save)

    df_mec = pd.DataFrame()
    df_result = pd.DataFrame()
    # 遍历基因组
    for genome in os.listdir(args.path_input_genome):
        print(genome)
        genome_name = genome.split(".")[0]
        # 建blast_db
        path_db_genome = make_blast_db_genome(genome, path_db_dir)
        # mec基因的比对，找到contig
        path_queries_mec = path_convert(os.path.join(args.scheme, "mec"))
        df_mec_alleles_genome = find_alleles(path_queries_mec, path_db_genome, path_blast_result, genome, 0.3)
        path_queries_ccr = path_convert(os.path.join(args.scheme, "ccr"))
        df_ccr_alleles_genome = find_alleles(path_queries_ccr, path_db_genome, path_blast_result, genome, 0.3)
        if df_mec.shape[0] == 0:
            df_mec = df_mec_alleles_genome
        else:
            df_mec = pd.concat([df_mec, df_mec_alleles_genome], axis=0)

        # 找出有mec的基因组
        if df_mec_alleles_genome.shape[0] == 1 and df_mec_alleles_genome.iloc[0]["contig"] == 0:
            # print(f"{genome}_no_mec")
            df_result_genome = pd.DataFrame({"genome": [genome_name], "scc_id": ["0"], "mec_complex": ["no_mec"], "note_complex":["no_mec"], "mec_coverage":["0"],
              "mecI_scc":["0"], "mecI_genome":["0"], "mecR_scc":["0"], "mecR_genome":["0"],
              "IS_left":["0"], "IS_right":["0"], "ccr_scc":["0"], "ccr_genome":["0"]})
            if df_result.shape[0] == 0:
                df_result = df_result_genome
            else:
                df_result = pd.concat([df_result, df_result_genome], axis=0)
        else:
            # 在全基因组中找mecI和mecR
            path_queries_meci = path_convert(os.path.join(args.scheme, "mecs"))
            df_meci_alleles_genome = find_alleles(path_queries_meci, path_db_genome, path_blast_result, genome, 0.01)

            # 根据mec位置找到scc，保存序列
            scc_path = extract_from_genome(genome, "scc_mec", df_mec_alleles_genome, 4000, 7000, path_scc_save)
            # 分析scc里的结构
            for scc_file in scc_path:
                # 建blast_db_scc
                path_db_scc = make_blast_db_scc(scc_file, path_db_dir)
                # mec基因的比对，找到contig
                path_queries_mec = path_convert(os.path.join(args.scheme, "mec"))
                df_mec_alleles_scc = find_alleles(path_queries_mec, path_db_scc, path_blast_result, scc_file, 0.3)
                df_meci_alleles_scc = find_alleles(path_queries_meci, path_db_scc, path_blast_result, scc_file, 0.01)
                path_queries_is = path_convert(os.path.join(args.scheme, "IS"))
                df_is_alleles_scc = find_alleles(path_queries_is, path_db_scc, path_blast_result, scc_file, 0.01)
                path_queries_ccr = path_convert(os.path.join(args.scheme, "ccr"))
                df_ccr_alleles_scc = find_alleles(path_queries_ccr, path_db_scc, path_blast_result, scc_file, 0.3)
                path_queries_ugpq = path_convert(os.path.join(args.scheme, "ugpQ"))
                df_ugpq_alleles_scc = find_alleles(path_queries_ugpq, path_db_scc, path_blast_result, scc_file, 0.3)

                extract_from_scc(scc_file, df_mec_alleles_scc, df_ugpq_alleles_scc, df_is_alleles_scc, path_mec_complex_save, path_hvr_save, path_except_hvr_save)
                # mec_complex_path, hvr_path =
                mec_id, type_mec_complex, note_complex = judge_mec_complex(df_mec_alleles_scc, df_is_alleles_scc, df_meci_alleles_scc, df_meci_alleles_genome)
                # judge_ccr_complex(df_ccr_alleles_scc, df_ccr_alleles_genome)
                df_result_genome = df_result_emerge(mec_id, type_mec_complex, note_complex, df_mec_alleles_scc, df_meci_alleles_scc, df_meci_alleles_genome, df_is_alleles_scc, df_ccr_alleles_scc, df_ccr_alleles_genome)
                if df_result.shape[0] == 0:
                    df_result = df_result_genome
                else:
                    df_result = pd.concat([df_result, df_result_genome], axis=0)

        # 删除genome的db文件
        if args.if_save_db == "1":
            execute(f"rm {path_db_genome}*")

    df_mec.to_excel("mec_alleles.xlsx", index=False)
    df_result.to_excel("mec_complex_typing_result.xlsx", index=False)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
进行MLST分型
使用之前确保当前环境已安装blast, python需求的包：executor, pandas, biopython
    ''', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-s', "--scheme", action="store", help="输入scheme的路径", default='./scheme')
    parser.add_argument("-g", "--path_input_genome", action="store", required=True,help="输入待比对基因组文件的路径")
    parser.add_argument("-o", "--output_dir", action="store", default=".", help="结果输出的路径")
    parser.add_argument("-db", "--if_save_db", action="store", default="1", help="是否保存db文件")
    parser.add_argument("-json", "--if_save_json", action="store", default="1", help="是否保存json文件")

    args = parser.parse_args()
    main()  # run main script