#!/bin/bash
echo "###################################################################################################"
echo " WELCOME TO CLUSTERING AND ANALAYSIS OF NEIGHBORHOOD GENES"
echo "PLEASE CHOOSE FROM THE FOLLOWING OPTIONS FOR FURTHER EXECUTION:"
echo "1)LOCAL .GBK AND .RGI ANNOTATED FILE EXISTS"
echo "2)ONLY LOCAL .FNA FILE EXISTS "
echo "###################################################################################################"



read -p "Please enter the option : " option

echo ${option}

case ${option} in 
     1)     read -p "Please enter the RGI file/files path : " RGI_Path
            read -p "Please enter the GBK file/files path : " GBK_Path


            cd $RGI_Path


            count_of_files=$(ls -1 | wc -l)
            echo $count_of_files


            cd "$OLDPWD"

            if [ $count_of_files -gt 1 ]
            then
              mkdir Multiple_instance
              mkdir Single_instance
              mkdir contigend_visualizations
              python git_upload_file_part1.py ${RGI_Path} ${GBK_Path}

              cd contigend_visualizations
              convert *.png -append complete_genome_visualizations.png


              cd Multiple_instance
              mkdir output_blast_multiple_instance
              for i in *.fasta
              do    
              
                     name=$(echo $i | cut -d'.' -f1)
                     echo $name
                     sed -i s/\'//g $i; #########removes single quotes from the file
                     sed -i 's/[][]//g' $i; ############removes square braces from the file

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
                     sed -i s/\'//g $i; #########removes single quotes from the file
                     sed -i 's/[][]//g' $i; ############removes square braces from the file

                     makeblastdb -dbtype prot -in $i
                     blastp -query $i -db $i -outfmt  "6 qseqid sseqid pident length evalue bitscore qseq sseq " -out output_blast_single_instance/$name.txt
                     
                    
              done

              cd "$OLDPWD"

              mkdir Multiple_instance_clusters
              mkdir Single_instance_clusters
              mkdir Phylogeny_fasta_files


              python git_upload_file_part2.py ${RGI_Path} ${GBK_Path} Single_instance/output_blast_single_instance Multiple_instance/output_blast_multiple_instance
          

            fi

            if [ $count_of_files -eq 1 ]
            then 
              mkdir contigend_visualizations_single_genome     
              python git_upload_file_part1.py ${RGI_Path} ${GBK_Path}                              
          
              cd contigend_visualizations_single_genome
              convert *.jpeg -append complete_genome_visualizations.jpeg

              cd "$OLDPWD"
            fi


          
            
            ;; 
   2) read -p "Please enter the .fna assembly files path : " Fna_Path_RGI_Annotation
      
      
      cd ${Fna_Path_RGI_Annotation}  

      mkdir allgbksrequired 
      mkdir allfnasrequired
      mkdir allrgisrequired


      CONDA_BASE=$($CONDA_EXE info --base)
      source $CONDA_BASE/etc/profile.d/conda.sh

      conda activate /home/chandana/anaconda3/envs/prokka/


      for i in *.fna
      do
         name=$(echo $i | cut -d'.' -f1)          
         prokka $i --prefix $name --outdir ${name}_prokka  --locustag $name    

         cp ${name}_prokka/${name}.gbk allgbksrequired
         cp ${name}_prokka/${name}.fna allfnasrequired
         
      done


      CONDA_BASE=$($CONDA_EXE info --base)
      source $CONDA_BASE/etc/profile.d/conda.sh

      conda activate /home/chandana/anaconda3/envs/rgi/   

      count_of_files=0
      for i in allfnasrequired/*.fna
      do
         name=$(basename $i | cut -d '.' -f1)
         echo $name         
         count_of_files=`expr $count_of_files + 1`
         rgi main --clean --input_sequence $i --alignment_tool BLAST  --num_threads 1 --output $name    
         
      done 
      
      for i in *.txt
      do
         
         mv $i allrgisrequired
         
      done


      echo $count_of_files

      cd "$OLDPWD" 
      ls

      if [ $count_of_files -gt 1 ]
      then
        mkdir Multiple_instance
        mkdir Single_instance
        mkdir contigend_visualizations


        python Neighborhood_Generator.py ${Fna_Path_RGI_Annotation}/allrgisrequired ${Fna_Path_RGI_Annotation}/allgbksrequired


        cd Multiple_instance
        mkdir output_blast_multiple_instance
        for i in *.fasta
        do     
        
               name=$(echo $i | cut -d'.' -f1)
               echo $name
               sed -i s/\'//g $i; #########removes single quotes from the file
               sed -i 's/[][]//g' $i; ############removes square braces from the file

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
               sed -i s/\'//g $i; #########removes single quotes from the file
               sed -i 's/[][]//g' $i; ############removes square braces from the file

               makeblastdb -dbtype prot -in $i
               blastp -query $i -db $i -outfmt  "6 qseqid sseqid pident length evalue bitscore qseq sseq " -out output_blast_single_instance/$name.txt
               
              
        done

        cd "$OLDPWD"

        mkdir Multiple_instance_clusters
        mkdir Single_instance_clusters
        mkdir Phylogeny_fasta_files


        
        python Cluster_Generator.py ${Fna_Path_RGI_Annotation}/allrgisrequired ${Fna_Path_RGI_Annotation}/allgbksrequired Single_instance/output_blast_single_instance Multiple_instance/output_blast_multiple_instance
       

      fi

      if [ $count_of_files -eq 1 ]
      then 
        mkdir contigend_visualizations_single_genome     
        python Neighborhood_Generator.py ${Fna_Path_RGI_Annotation}/allrgisrequired ${Fna_Path_RGI_Annotation}/allgbksrequired
                            
    
        cd contigend_visualizations_single_genome
        convert *.png -append complete_genome_visualizations.png

        cd "$OLDPWD"
      fi

      echo "###################################################################################################"
      echo "THANK YOU"
      exit 1 
      ;;
      
  
esac 
