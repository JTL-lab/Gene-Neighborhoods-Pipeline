

echo " welcome to clustering and analaysis of neighborhood genes"
echo "Please choose from the following options for further execution:"
echo "1)Local .gbk and .rgi annotated file exists"
echo "2)Only local .fna file exists "


read -p "Please enter the option : " option

echo ${option}

case ${option} in 
     1)     read -p "Please enter the RGI file/files path : " RGI_Path
            read -p "Please enter the GBK file/files path : " GBK_Path
            
            python git_upload_file_part1.py ${RGI_Path} ${GBK_Path}

            mkdir Multiple_instance
            mkdir Single_instance

            mkdir contigend_visualizations
            
            


            
            for i in *.png
            do
               mv $i contigend_visualizations/$i
            done



           



            ### creating separate directory for .fasta files obtain from python script and running ALL-Vs-All Blast to get similarity results.
            
            cd Multiple_instance
            mkdir output_blast_multiple_instance
            for i in *.fasta
            do    
            
                   name=$(echo $i | cut -d'.' -f1)
                   echo $name
                   makeblastdb -dbtype prot -in $i
                   blastp -query $i -db $i -outfmt  "6 qseqid sseqid pident length evalue bitscore qseq sseq " -out output_blast_multiple_instance/$name.txt
                   #mv $i.txt output_blast/$name.txt
                  
            done
            
            cd "$OLDPWD"

            cd Single_instance
            mkdir output_blast_single_instance

            for i in *.fasta
            do    
            
                   name=$(echo $i | cut -d'.' -f1)
                   echo $name
                   makeblastdb -dbtype prot -in $i
                   blastp -query $i -db $i -outfmt  "6 qseqid sseqid pident length evalue bitscore qseq sseq " -out output_blast_single_instance/$name.txt
                   #mv $i.txt output_blast/$name.txt
                  
            done

            cd "$OLDPWD"

            mkdir Multiple_instance_clusters
            mkdir Single_instance_clusters
            mkdir Phylogeny_fasta_files


            python git_upload_file_part2.py Fin_RGI Fin_gbk Single_instance/output_blast_single_instance Multiple_instance/output_blast_multiple_instance
          
            


            ;; 
   2) read -p "Please enter the .fna assembly files path : " Fna_Path_RGI_Annotation
      
      
      # cd ${Fna_Path_RGI_Annotation}     

      # mkdir allgbksrequired

      # mkdir allrgisrequired
 

      # # for i in *.fna
      # # do
      # #    name=$(echo $i | cut -d'.' -f1) 
      # #    echo $i

      # #    prokka $i --prefix $name --outdir ${name}_prokka  --locustag $name


      # #    rgi main --clean --input_sequence $i --alignment_tool BLAST  --num_threads 1 --output $name

      # #    cp ${name}_prokka/${name}.gbk allgbksrequired
         
      # # done


      # for i in *.fna
      # do
      #    name=$(echo $i | cut -d'.' -f1) 
      #    echo $i      


      #    rgi main --clean --input_sequence $i --alignment_tool BLAST  --num_threads 1 --output $name

         
      #    cp ${name}.txt allrgisrequired
      # done

      # cd "$OLDPWD"



      python git_upload_file_part1.py ${Fna_Path_RGI_Annotation}/allrgisrequired ${Fna_Path_RGI_Annotation}/allgbksrequired 



      mkdir Multiple_instance
      mkdir Single_instance

      mkdir contigend_visualizations     
      


      
      for i in *.png
      do
         mv $i contigend_visualizations/$i
      done



     



      ### creating separate directory for .fasta files obtain from python script and running ALL-Vs-All Blast to get similarity results.
      
      cd Multiple_instance
      mkdir output_blast_multiple_instance
      for i in *.fasta
      do    
      
             name=$(echo $i | cut -d'.' -f1)
             echo $name
             makeblastdb -dbtype prot -in $i
             blastp -query $i -db $i -outfmt  "6 qseqid sseqid pident length evalue bitscore qseq sseq " -out output_blast_multiple_instance/$name.txt
             #mv $i.txt output_blast/$name.txt
            
      done
      
      cd "$OLDPWD"

      cd Single_instance
      mkdir output_blast_single_instance

      for i in *.fasta
      do    
      
             name=$(echo $i | cut -d'.' -f1)
             echo $name
             makeblastdb -dbtype prot -in $i
             blastp -query $i -db $i -outfmt  "6 qseqid sseqid pident length evalue bitscore qseq sseq " -out output_blast_single_instance/$name.txt
             #mv $i.txt output_blast/$name.txt
            
      done

      cd "$OLDPWD"

      mkdir Multiple_instance_clusters
      mkdir Single_instance_clusters
      mkdir Phylogeny_fasta_files


      python git_upload_file_part2.py ${Fna_Path_RGI_Annotation}/allrgisrequired ${Fna_Path_RGI_Annotation}/allgbksrequired Single_instance/output_blast_single_instance Multiple_instance/output_blast_multiple_instance
    

    
      exit 1 
      ;;
      
  
esac 
