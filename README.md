# asexp
Identify and calculate  the expression of alternative splicing events
各个脚本的具体使用可以通过 'perl 脚本名称 -h' 查看帮助信息
1.Treat_mapping.pl: 对比对后的bam文件筛选，输出reads的跨越的junction信息，missmatch信息及overhang信息等
2.Filter_gtf.pl: 通过第一步输出的信息，对组装出来的gtf文件中转录本进行筛选，获得高可信度的转录本信息
3.Discovery_AS.pl: 对筛选后的gtf格式的转录本文件，进行可变剪接事件的识别和分类
4.AS_exp.pl: 通过第一步生成的reads junction信息等，对第三步识别的可变剪接事件定量
