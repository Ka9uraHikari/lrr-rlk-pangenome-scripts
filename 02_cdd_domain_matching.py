import csv


cddid_file = "/path/to/data/CDDdatabase/cddid_all.tbl"
scan_file = "/path/to/project/07_CDDscan/all_cdd_named.out"
output_file = "/path/to/project/07_CDDscan/all_cdd_named_full.out"


cdd_map = {}
with open(cddid_file, newline='') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        if len(row) < 4:
            continue
        cdd_id = row[0].strip()
        cdd_map[cdd_id] = [row[1].strip(), row[2].strip(), row[3].strip()]


with open(scan_file, newline='') as f_in, open(output_file, 'w', newline='') as f_out:
    reader = csv.reader(f_in, delimiter='\t')
    writer = csv.writer(f_out, delimiter='\t')

    for row in reader:
        if not row:
            continue
        if row[0].startswith("#"):
            writer.writerow(row)
            continue
        cdd_field = row[1].strip()
        cdd_id = cdd_field.replace("CDD:", "")
        if cdd_id in cdd_map:
            row.extend(cdd_map[cdd_id])
        writer.writerow(row)

print(f"完成：生成带可读ID的文件 {output_file}")
