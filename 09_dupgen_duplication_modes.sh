#!/bin/bash
set -e

GENOME_DIR=/path/to/project/annotation
OUT_DIR=/path/to/project/14_复制类型/dupgen
ARAB=/path/to/project/14_复制类型/Araport11_pep_20250411.fa
SCRIPT=/path/to/data/dupgen/DupGen_finder-1.0.0/DupGen_finder.pl

THREADS=20
mkdir -p "$OUT_DIR"

awk 'BEGIN{OFS=""}
/^>/{
    id=$0
    sub(/^>/,"",id)
    sub(/ .*/,"",id)
    print ">" id
    next
}
{print}' "$ARAB" > "$OUT_DIR/Arabidopsis.pep"

awk 'BEGIN{OFS="\t"; n=0}
    /^>/{
        id=$0
        sub(/^>/,"",id)
        n++
        print "chrAT", id, n*1000, n*1000+999
    }
' "$OUT_DIR/Arabidopsis.pep" > "$OUT_DIR/Arabidopsis.gff"

for pep in $GENOME_DIR/*_pep.fa; do
    base=$(basename "$pep")
    name=${base%%_pep.fa}

    echo "========================"
    echo "Processing $name"
    echo "========================"

    gff=$GENOME_DIR/${name}.gff3
    if [ ! -f "$gff" ]; then
        echo "❌ Missing GFF for $name"
        continue
    fi

    mkdir -p "$OUT_DIR/$name"
    cd "$OUT_DIR/$name"

    awk 'BEGIN{OFS=""}
    /^>/{
        id=$0
        sub(/^>/,"",id)
        sub(/ .*/,"",id)
        sub(/\.[0-9]+$/,"",id)
        print ">" id
        next
    }
    {print}' "$pep" > "${name}.pep"

    awk 'BEGIN{OFS="\t"}
        $3=="gene"{
            attr=$9
            if (match(attr,/ID=[^;]+/)) {
                id = substr(attr, RSTART+3, RLENGTH-3)
                sub(/\.[0-9]+$/,"",id)
                print $1, id, $4, $5
            }
        }' "$gff" > "${name}.gff"

    diamond makedb --in "${name}.pep" -d "${name}"
    diamond blastp -d "${name}" -q "${name}.pep" -o "${name}.blast" \
        -e 1e-5 -k 5 -p "$THREADS" --quiet

    diamond makedb --in "$OUT_DIR/Arabidopsis.pep" -d Arabidopsis
    diamond blastp -d Arabidopsis -q "${name}.pep" -o "${name}_Arabidopsis.blast" \
        -e 1e-5 -k 5 -p "$THREADS" --quiet

    cat "${name}.gff" "$OUT_DIR/Arabidopsis.gff" > "${name}_Arabidopsis.gff"

    perl "$SCRIPT" \
        -i . \
        -t "${name}" \
        -c Arabidopsis \
        -o result

done
