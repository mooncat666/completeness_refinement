require(dplyr)
require(stringr)

gaf=commandArgs(trailingOnly=TRUE)[1]
par=commandArgs(trailingOnly=TRUE)[2]

df <- read.csv(gaf,sep='\t', header=FALSE)
colnames(df) <- c('contig', 'contig_length', 'q_start','q_end','q_strand','path','sth1','r_start','r_end','sth2',
                  'alignment_length', 'mapq','mismatches','ali_score','divergence','identity', 'cigar')
df <- mutate(df, identity=as.numeric(as.character(str_replace(identity,"id:f:",""))))
df[c('contig_name','flag','coverage','len')] <- str_split_fixed(df$contig, ' ', 4)
df <- mutate(df,coverage=as.numeric(as.character(str_replace(coverage,"multi=",""))))

count_cov_length <- function(df) {
    length=0
    x=df[1,'q_start']
    y=df[1,'q_end']
    if (nrow(df)>1) {
        for (row in 2:nrow(df)) {
            if (df[row,'q_start']>y) { # no overlap
                # add length
                length=length+(y-x)
                x=df[row,'q_start']
                y=df[row,'q_end']
            } else { # overlap
                if (df[row,'q_end']>y) { # not a containment
                    y=df[row,'q_end']
                }
            }
        }
    }
    length=length+(y-x)
    return(length)
}


sm <- group_by(df,contig_name) %>% 
    filter(alignment_length>50) %>%
    arrange(q_start,.by_group=TRUE) %>%
    summarise(cov_length=count_cov_length(cur_data()),contig_length=mean(contig_length),coverage=mean(coverage), identity=mean(identity)) %>%
    mutate(perc_covered=cov_length/contig_length) %>%
    filter(perc_covered>par)


names <- select(sm,contig_name)


write.table(names,stdout(),row.names=FALSE,col.names=FALSE)
