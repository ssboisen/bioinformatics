(ns rosalind.core
  (:require [clojure.string :as string]))

(defn output-to-file [file line-vector]
  (spit file
        (apply str (interpose "\n" line-vector))))

(defn parse-fasta [input]
  (let [fasta-parts (rest (string/split input #">"))
        pairs (for [part fasta-parts
                    :let [[name & dnastring] (string/split part #"\n")]]
                [name (string/join "" dnastring)])]
    pairs))

(defn partition-dna [dnastring]
  (map (partial apply str) (partition 3 dnastring)))

(defn is-sorted?
  ([coll]
   (is-sorted? compare coll))
  ([^java.util.Comparator comp coll]
   (every? (partial apply comp) (partition 2 1 coll))))

(def monoisotopic-mass-table
  {
   "A"   71.03711
   "C"   103.00919
   "D"   115.02694
   "E"   129.04259
   "F"   147.06841
   "G"   57.02146
   "H"   137.05891
   "I"   113.08406
   "K"   128.09496
   "L"   113.08406
   "M"   131.04049
   "N"   114.04293
   "P"   97.05276
   "Q"   128.05858
   "R"   156.10111
   "S"   87.03203
   "T"   101.04768
   "V"   99.06841
   "W"   186.07931
   "Y"   163.06333
   })

(def dna-codon-table
    {
    "TTT" "F"      "CTT" "L"      "ATT" "I"      "GTT" "V"
    "TTC" "F"      "CTC" "L"      "ATC" "I"      "GTC" "V"
    "TTA" "L"      "CTA" "L"      "ATA" "I"      "GTA" "V"
    "TTG" "L"      "CTG" "L"      "ATG" "M"      "GTG" "V"
    "TCT" "S"      "CCT" "P"      "ACT" "T"      "GCT" "A"
    "TCC" "S"      "CCC" "P"      "ACC" "T"      "GCC" "A"
    "TCA" "S"      "CCA" "P"      "ACA" "T"      "GCA" "A"
    "TCG" "S"      "CCG" "P"      "ACG" "T"      "GCG" "A"
    "TAT" "Y"      "CAT" "H"      "AAT" "N"      "GAT" "D"
    "TAC" "Y"      "CAC" "H"      "AAC" "N"      "GAC" "D"
    "TAA" "Stop"   "CAA" "Q"      "AAA" "K"      "GAA" "E"
    "TAG" "Stop"   "CAG" "Q"      "AAG" "K"      "GAG" "E"
    "TGT" "C"      "CGT" "R"      "AGT" "S"      "GGT" "G"
    "TGC" "C"      "CGC" "R"      "AGC" "S"      "GGC" "G"
    "TGA" "Stop"   "CGA" "R"      "AGA" "R"      "GGA" "G"
    "TGG" "W"      "CGG" "R"      "AGG" "R"      "GGG" "G"
     })

(def rna-codon-table
  {
     "UUU" "F"      "CUU" "L"      "AUU" "I"      "GUU" "V"
     "UUC" "F"      "CUC" "L"      "AUC" "I"      "GUC" "V"
     "UUA" "L"      "CUA" "L"      "AUA" "I"      "GUA" "V"
     "UUG" "L"      "CUG" "L"      "AUG" "M"      "GUG" "V"
     "UCU" "S"      "CCU" "P"      "ACU" "T"      "GCU" "A"
     "UCC" "S"      "CCC" "P"      "ACC" "T"      "GCC" "A"
     "UCA" "S"      "CCA" "P"      "ACA" "T"      "GCA" "A"
     "UCG" "S"      "CCG" "P"      "ACG" "T"      "GCG" "A"
     "UAU" "Y"      "CAU" "H"      "AAU" "N"      "GAU" "D"
     "UAC" "Y"      "CAC" "H"      "AAC" "N"      "GAC" "D"
     "UAA" "Stop"   "CAA" "Q"      "AAA" "K"      "GAA" "E"
     "UAG" "Stop"   "CAG" "Q"      "AAG" "K"      "GAG" "E"
     "UGU" "C"      "CGU" "R"      "AGU" "S"      "GGU" "G"
     "UGC" "C"      "CGC" "R"      "AGC" "S"      "GGC" "G"
     "UGA" "Stop"   "CGA" "R"      "AGA" "R"      "GGA" "G"
     "UGG" "W"      "CGG" "R"      "AGG" "R"      "GGG" "G"
   })

(def molecule-weight-table
  {
   "G" 57
   "A" 71
   "S" 87
   "P" 97
   "V" 99
   "T" 101
   "C" 103
   "I" 113
   "L" 113
   "N" 114
   "D" 115
   "K" 128
   "Q" 128
   "E" 129
   "M" 131
   "H" 137
   "F" 147
   "R" 156
   "Y" 163
   "W" 186
   })

(def reverse-molecule-weight-table (into {} (map (fn [[k v]] [v k]) molecule-weight-table)))
(def dna-codon-to-protein (partial get dna-codon-table))

(def rna-codon-to-protein (partial get rna-codon-table))

(defn partition-string
  ([n s]
     (partition-string n n s))
  ([n step s]
     (loop [start 0
            acc []]
       (if (> (+ start n) (count s))
         acc
         (recur (+ start step) (conj acc (subs s start (+ start n))))))))

(defn abs [x]
  (if (< x 0)
    (* x -1)
    x))
