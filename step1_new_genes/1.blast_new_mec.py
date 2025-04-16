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


def main():
    # 建立临时文件夹和输出文件夹
    # 结果文件夹
    path_result_dir = path_convert(os.path.join(args.output_dir, args.name_output_dir))
    exist_mk(path_result_dir)

    # coverage低于阈值的保存位置
    path_new_result = path_convert(os.path.join(path_result_dir, "3.new_genes"))
    exist_mk(path_new_result)
    # blast结果保存
    path_blast_result = path_convert(os.path.join(path_result_dir, "2.blast_result"))
    exist_mk(path_blast_result)
    # blast_db保存
    path_db_dir = path_convert(os.path.join(path_result_dir, "1.db"))
    exist_mk(path_db_dir)

    result_all = pd.DataFrame()
    # 遍历基因组
    for genome in os.listdir(args.path_input_genome):
        # 建blast_db
        genome_name = genome.split(".")[0]
        title = f"blastdb_{genome_name}"
        path_db = path_convert(os.path.join(path_db_dir, title))
        path_title = path_convert(os.path.join(path_db_dir, f"blastdb_{genome_name}"))
        path_genome = path_convert(os.path.join(args.path_input_genome, genome))
        make_blast_db(path_genome, path_title, path_db)

        result_df = pd.DataFrame()

        mec_genes = ["mecA", "mecB", "mecC", "mecD", "mecA1", "mecA2", "mecC1", "mecC2"]

        for gene_name in mec_genes:
            # 进行blast
            blast_n(gene_name, genome_name, path_blast_result, path_title)
            # 解析blast结果
            df_hsp = parse_blast_json(gene_name, genome_name, path_blast_result)
            # 判断有没有结果，有的话继续
            if not df_hsp.empty:
                # 判断是否为第一个有内容的结果，是的话直接存结果
                if result_df.empty:
                    result_df = df_hsp
                else:
                    # 结果中进行比较，先根据position确认比对为同一个位置，然后比较，更新
                    for index, row_a in df_hsp.iterrows():
                        df_contig = result_df[result_df.pos_contig == row_a["pos_contig"]]
                        if df_contig.empty:
                            # 增加这一行
                            result_df.loc[len(result_df)] = row_a
                            # print("new_row")
                        else:
                            # 依次判断是否为同一position
                            contain_row = 0
                            for index_b, row_b in df_contig.iterrows():
                                if judge_pos(row_a, row_b):
                                    contain_row = 1
                                    if float(row_a["coverage"]) > float(row_b["coverage"]):
                                        result_df.loc[index_b] = row_a
                            # 判断没有同一position的row
                            if contain_row == 0:
                                result_df.loc[len(result_df)] = row_a
        if result_df.empty:
            result_output = {"genome": [genome_name], "query_name": ["empty"], "hit_count": [0], "hsp_count": [0],
                             "coverage": [0.0], "identity": [0.0], "pos_contig": [0], "pos_from": [0], "pos_to": [0],
                             "strand": ["Plus"], "hseq": ["A"]}
            result_df = pd.DataFrame(result_output)
        # 删除db文件
        if args.if_save_db == "1":
            execute(f"rm {path_title}*")
        # 获得每个基因组的result_df
        result_all = pd.concat([result_all, result_df], axis=0)
    print(result_all)
    for index, row in result_all.iterrows():
        if 0.1 < float(row["coverage"]) < float(args.coverage_threshold) and float(row["identity"]) < 0.95:
            save_fa(row, path_new_result)
        pass

    file_saved_all = path_convert(os.path.join(path_result_dir, f"result_all_{args.name_output_dir}.tsv"))
    result_all.to_csv(file_saved_all, sep="\t", index=False)


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
    # if os.path.isdir(p):
    #     if os.path.exists(p):
    #         pass
    #         print(2)
    #     else:
    #         print(3)
    #         makedir(p)


def makedir(p):
    if os.path.exists(p):
        # print("文件已存在")
        pass
    else:
        execute(f"mkdir {p}")


def make_blast_db(input_genome, title, path_output):
    # 建立blast数据库
    execute(f"makeblastdb -in {input_genome} -dbtype nucl -title {title} -parse_seqids -out {path_output}")

def blast_n(gene_name, genome_name, path_blast_result, path_db):
    # 进行blast操作
    path_gene = path_convert(os.path.join(args.query_gene, f"{gene_name}.fa"))
    path_output = path_convert(os.path.join(path_blast_result, f"{genome_name}_{gene_name}.json"))

    execute(f"blastn -query {path_gene} -out {path_output} -outfmt 15 -db {path_db} -evalue 0.001")

def parse_blast_json(gene_name, genome_name, path_blast_result):
    json_path = path_convert(os.path.join(path_blast_result, f"{genome_name}_{gene_name}.json"))
    file = open(json_path)
    f = json.load(file)
    file.close()

    # 构建primary result
    result_output = dict()
    result_output["genome"] = []
    result_output["query_name"] = []
    result_output["hit_count"] = []
    result_output["hsp_count"] = []
    result_output["coverage"] = []
    result_output["identity"] = []
    result_output["pos_contig"] = []
    result_output["pos_from"] = []
    result_output["pos_to"] = []
    result_output["strand"] = []
    result_output["hseq"] = []
    querys = f["BlastOutput2"]
    for query in querys:
        results = query["report"]["results"]["search"]
        query_len = results["query_len"]
        hits = results["hits"]
        hit_num = 0
        for hit in hits:
            hit_id = hit["description"][0]["id"]
            sum_identity = 0
            hsp_count = 0
            for hsp in hit["hsps"]:
                hsp_count += 1
                sum_identity += hsp["identity"]
            coverage = float("{:.2f}".format(sum_identity / query_len))
            hsp_max = hit["hsps"][0]
            hsp_hseq = hsp_max["hseq"]
            hsp_max_length = abs(hsp_max["hit_from"] - hsp_max["hit_to"])
            identity_max_hsp = float("{:.2f}".format((hsp_max_length - hsp_max["gaps"])/ hsp_max_length))

            if coverage > 0.2:
                hit_num += 1
                result_output["coverage"].append(coverage)
                result_output["identity"].append(identity_max_hsp)
                result_output["hit_count"].append(hit_num)
                result_output["hsp_count"].append(hsp_count)
                result_output["pos_from"].append(hsp_max["hit_from"])
                result_output["pos_to"].append(hsp_max["hit_to"])
                result_output["strand"].append(hsp_max["hit_strand"])
                result_output["genome"].append(hit_id.split(".")[0])
                result_output["pos_contig"].append(hit_id.split("ctg")[1])
                result_output["query_name"].append(gene_name)
                result_output["hseq"].append(hsp_hseq)
    df_result = pd.DataFrame(result_output)
    return df_result


def judge_pos(row1, row2):
    from_1 = min(int(row1["pos_from"]), int(row1["pos_to"]))
    to_1 = max(int(row1["pos_from"]), int(row1["pos_to"]))
    from_2 = min(int(row2["pos_from"]), int(row2["pos_to"]))
    to_2 = max(int(row2["pos_from"]), int(row2["pos_to"]))
    if abs(from_1 - from_2) < 100 and abs(to_1 - to_2) < 100:
        return True
    else:
        return False


def save_fa(row, path_new_result):
    query_name = row["query_name"]
    genome_name = row["genome"]
    hit_count = row["hit_count"]
    coverage = row["coverage"]
    hsp_count =row["hsp_count"]
    file_name = f"{query_name}_{genome_name}_hit{hit_count}_hsp{hsp_count}_{coverage}"
    pos_contig = row["pos_contig"]
    pos_from = row["pos_from"]
    pos_to = row["pos_to"]
    des = f"contig_{pos_contig}_from_{pos_from}_to_{pos_to}"
    file_path = path_convert(os.path.join(path_new_result, file_name))
    rec1 = SeqRecord(seq=Seq(row["hseq"]), id=file_name, description=des)
    myrecords = [rec1]
    SeqIO.write(myrecords, file_path, "fasta")



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
进行MLST分型
使用之前确保当前环境已安装blast, python需求的包：executor, pandas, biopython
    ''', formatter_class=RawTextHelpFormatter)
    parser.add_argument('-q', "--query_gene", action="store", required=True, help="输入scheme的路径")
    parser.add_argument("-g", "--path_input_genome", action="store", required=True,help="输入待比对基因组文件的路径")
    parser.add_argument("-o", "--output_dir", action="store", default=".", help="结果输出的路径")
    parser.add_argument("-on", "--name_output_dir", action="store", default=".", help="结果输出的路径文件夹名")
    parser.add_argument("-db", "--if_save_db", action="store", default="1", help="是否保存db文件")
    parser.add_argument("-json", "--if_save_json", action="store", default="1", help="是否保存json文件")
    parser.add_argument("-t", "--coverage_threshold", action="store", help="coverage判断标准")

    args = parser.parse_args()
    main()  # run main script