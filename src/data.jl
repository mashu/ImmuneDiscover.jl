module data
    using FASTX
    """
    Save records with sequences and quality scores to fastq file
    """
    function write_fastq(path, records)
        FASTQ.Writer(open(path, "w")) do writer
            for record in records
                write(writer, record)
            end
        end
    end
end