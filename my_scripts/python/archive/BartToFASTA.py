from Bio import SeqIO
import csv
import pandas as pd

def extract_transcripts(fasta_file, transcript_list, output_file):
    # Чтение всех транскриптов из FASTA
    transcripts = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    
    # Создание списка последовательностей, которые нужно извлечь
    selected_transcripts = [transcripts[transcript] for transcript in transcript_list if transcript in transcripts]
    
    # Запись выбранных транскриптов в новый FASTA файл
    SeqIO.write(selected_transcripts, output_file, "fasta")
    print(f"Extracted {len(selected_transcripts)} transcripts to {output_file}")

# Пример использования
fasta_file = "/Users/ПК/Documents/bachelor/data/reference_data/BaRT2v18.fasta" # Путь к исходному FASTA файлу
transcript_csv = pd.read_csv("/Users/ПК/Documents/bachelor/data/4_DataAnalysis/kallisto/UpRegGamma.csv") # Список транскриптов, которые нужно извлечь
transcript_list = transcript_csv.T.iloc[1].to_list()
output_file = "/Users/ПК/Documents/bachelor/data/4_DataAnalysis/kallisto/UpRegGamma.fasta"  # Путь для сохранения результатов

extract_transcripts(fasta_file, transcript_list, output_file)
