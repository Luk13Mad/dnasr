//cargo run generate --db false -i test.fasta -o test
extern crate bio;
use clap::ArgMatches;
use search::search_main;
use shannon_entropy::entropy_calculation_main;
#[macro_use]
extern crate itertools;

//structures

//Generate_Config struct
#[derive(Clone)]
pub struct Generate_Config {
    pub input_file: String,
    /*pub meth: bool,*/
    pub output_dir: String,
    pub window: u32,
    pub kmer: u32,
    pub np: u32,
    pub matrix_dim: u32,
    pub sub_matrix: u32,
    pub database: bool,
    pub offset: u32,
    pub early_exit: bool,
}
//Generate_Config implementation block
impl Generate_Config {
    pub fn new(args: &ArgMatches) -> Result<Generate_Config, &'static str> {
        let input_file = String::from(args.value_of("input_file").unwrap());
        /*let meth = args.value_of("meth").unwrap().parse::<bool>().unwrap();*/
        let output_dir = String::from(args.value_of("output_dir").unwrap());
        let window = args.value_of("window").unwrap().parse::<u32>().unwrap();
        let kmer = args.value_of("kmer").unwrap().parse::<u32>().unwrap() + 1;
        let np = args.value_of("np").unwrap().parse::<u32>().unwrap();
        let matrix_dim = args.value_of("matrix_dim").unwrap().parse::<u32>().unwrap();
        let sub_matrix = args.value_of("sub_matrix").unwrap().parse::<u32>().unwrap();
        let database = args.value_of("database").unwrap().parse::<bool>().unwrap();
        let offset = args.value_of("offset").unwrap().parse::<u32>().unwrap();
        let early_exit = args
            .value_of("early_exit")
            .unwrap()
            .parse::<bool>()
            .unwrap();
        let myconfig = Generate_Config {
            input_file,
            output_dir,
            window,
            kmer,
            np,
            matrix_dim,
            sub_matrix,
            database,
            offset,
            early_exit,
        };
        match myconfig.isValid() {
            Ok(()) => return Ok(myconfig),
            Err(e) => return Err(e),
        }
    }

    fn isValid(&self) -> Result<(), &'static str> {
        if self.window <= 0 || self.kmer >= self.window {
            return Err("window must be positive,window must be larger than kmer");
        } else if self.kmer < 3 {
            return Err("kmer must be larger or equal 3");
        } else if self.np <= 0 {
            return Err("at least one core must be specified");
        } else if self.matrix_dim <= 0 || self.sub_matrix >= self.matrix_dim {
            return Err("matrix_dim must be positive, sub_matrix must be smaller than matrix_dim");
        } else if self.sub_matrix < 2 {
            return Err("sub_matrix must be at least two");
        } else if self.offset < 1 {
            return Err("offset must be at least 1");
        } else if !(self.input_file.contains("fasta") | self.input_file.contains("fastq")) {
            return Err("file must be *.fasta or *.fastq");
        } else if Generate_Config::check_dimensions(self.kmer,self.matrix_dim,self.sub_matrix) {
            return Err("Dimensions sub_matrix, matrix_dim or kmer incompatible.")
        } else {
            return Ok(());
        }
    }

    fn check_dimensions(kmer : u32, mdim: u32, sub : u32)-> bool{
        let base : u32 = 4;
        if (mdim % sub == 0){
            for k in (3..kmer){
                if (base.pow(k) % sub != 0){
                    return true;
                }
            }
            return false;
        } else {
            return true;
        }
    }
}

//struct with sliding windows and associated entropy vectors
//struct

//Search_Config struct
#[derive(Clone)]
pub struct Search_Config {
    pub query: String,
    pub database: String,
    pub output_dir: String,
}
//Search_Config implementation block
impl Search_Config {
    pub fn new(args: &ArgMatches) -> Result<Search_Config, &'static str> {
        let query = String::from(args.value_of("query").unwrap());
        let database = String::from(args.value_of("database").unwrap());
        let output_dir = String::from(args.value_of("output_dir").unwrap());
        let myconfig = Search_Config {
            query,
            database,
            output_dir,
        };
        return Ok(myconfig);
    }
}

//functions
pub fn run_generate(config: Result<Generate_Config, &'static str>) -> Result<(), &'static str> {
    println!("generate will be run");
    match config {
        Ok(c) => entropy_calculation_main(c),
        Err(e) => return Err(e),
    }
}

pub fn run_search(config: Result<Search_Config, &'static str>) -> Result<(), &'static str> {
    println!("search will be run");
    match config {
        Ok(c) => search_main(c),
        Err(e) => return Err(e),
    }
}

//modules

//module calculate entropy
pub mod shannon_entropy {
    use super::*;
    use bio::io::fasta;
    use bio::io::fastq;
    use ndarray::{Array2,ShapeBuilder,Array,s};
    use rayon::prelude::*;
    use itertools::Itertools;
    use std::{
        collections::HashMap,
        fs::File,
        io::{BufWriter, Write},
        path::Path,
        process, str,
        vec::Vec,
    };

    //structs inside shannon_entropy module
    #[derive(Clone, Debug)]
    pub struct Entropy {
        pub kmers: Vec<u32>,
        pub kmer_dict: HashMap<String, Vec<String>>,
        pub entropy_values: HashMap<String, Vec<f64>>,
    }
    impl Entropy {
        fn new(_window: &str, kmers: &Vec<u32>) -> Result<Entropy, &'static str> {
            //let bases=String::from("ATGC");
            let kmers = kmers.clone();
            let mut kmer_dict = HashMap::new();
            for k in &kmers {
                kmer_dict.insert(k.to_string(), Entropy::kproduct(String::from("ATGC"), *k));
            }
            let mut entropy_map = HashMap::new();
            for k in kmer_dict.values() {
                for l in k {
                    entropy_map.insert(l.clone(), vec![]);
                }
            }
            return Ok(Entropy {
                kmers,
                kmer_dict,
                entropy_values: entropy_map,
            });
        }

        fn kproduct(seq: String, k: u32) -> Vec<String> {
            match k {
                0 => vec![],
                1 => seq.chars().map(|c| c.to_string()).collect(),
                2 => iproduct!(seq.chars(), seq.chars())
                    .map(|(a, b)| format!("{}{}", a, b))
                    .collect(),
                _ => iproduct!(Entropy::kproduct(seq.clone(), k - 1), seq.chars())
                    .map(|(a, b)| format!("{}{}", a, b))
                    .collect(),
            }
        }
    }

    #[derive(Clone, Debug)]
    pub struct Record {
        pub seq: String,
        pub id: String,
        /*pub meth : bool,
        pub meth_scale : f32,*/
        pub window_size: u32,
        pub window_starts: Vec<String>,
        pub kmers: Vec<u32>,
        pub entropy_vectors: Vec<Entropy>,
    }
    impl Record {
        fn new(config: &Generate_Config, seq: &[u8], id: &str) -> Result<Record, &'static str> {
            let seq = std::str::from_utf8(seq.clone()).unwrap().to_string();
            let id = String::from(id);
            /*let meth = config.meth;
            let meth_scale=match meth{true => 1.0,false => 1.0,}; //TODO implement if meth == true*/
            let window_size = config.window;
            let kmers = (3..config.kmer).collect();
            let window_starts = match Record::generate_window_starts(&seq, &window_size) {
                Ok(r) => r,
                Err(_e) => return Err("Error splitting seq in windows"),
            };
            let entropy_vectors = match Record::generate_entropy_vectors(&kmers, &window_starts) {
                Ok(r) => r,
                Err(_e) => return Err("Error generating entropy vectors"),
            };
            return Ok(Record {
                seq,
                id,
                window_size,
                window_starts,
                kmers,
                entropy_vectors,
            });
        }

        fn generate_window_starts<'a, 'b>(
            seq: &'a str,
            window_size: &'b u32,
        ) -> Result<Vec<String>, &'static str> {
            let window_vector_bytes: Vec<&[u8]> =
                seq.as_bytes().windows(*window_size as usize).collect();
            let window_vector_string: Vec<String> = window_vector_bytes
                .iter()
                .map(|x| match str::from_utf8(x) {
                    Ok(v) => String::from(v),
                    Err(_e) => panic!("Invalid UTF-8 sequence in input file"),
                })
                .collect();
            return Ok(window_vector_string);
        }

        fn generate_entropy_vectors(
            kmers: &Vec<u32>,
            window_starts: &Vec<String>,
        ) -> Result<Vec<Entropy>, &'static str> {
            let mut out_vector = Vec::new();
            for window in window_starts {
                out_vector.push(match Entropy::new(window, kmers) {
                    Ok(v) => v,
                    Err(_e) => return Err("Error in Entropy::new()"),
                });
            }
            return Ok(out_vector);
        }

        fn get_position(&mut self) -> () {
            for (a, b) in self
                .window_starts
                .iter()
                .zip(self.entropy_vectors.iter_mut())
            {
                for s in &self.kmers {
                    let mychunks = Record::char_windows(a, *s as usize).collect::<Vec<&str>>();
                    let myindex = (0..mychunks.len() as u32).collect::<Vec<u32>>();
                    for (chunk, idx) in mychunks.iter().zip(myindex.iter()) {
                        let position_vec = match b.entropy_values.get_mut(*chunk) {
                            Some(r) => r,
                            None => continue,
                        };
                        position_vec.push(*idx as f64);
                    }
                }
            }
            return ();
        }

        fn char_windows<'a>(src: &'a str, win_size: usize) -> impl Iterator<Item = &'a str> {
            src.char_indices().flat_map(move |(from, _)| {
                src[from..]
                    .char_indices()
                    .skip(win_size - 1)
                    .next()
                    .map(|(to, c)| &src[from..from + to + c.len_utf8()])
            })
        }

        fn calc_alpha(&mut self) -> () {
            for a in self.entropy_vectors.iter_mut() {
                for (_, val) in a.entropy_values.iter_mut() {
                    if val.len() != 0 {
                        let mut first_el: Vec<f64> = vec![1.0 / val[0]];
                        let mut rest: Vec<f64> = Vec::new();
                        for (x, y) in val[1..].iter().zip(val.iter()) {
                            rest.push(1.0 / (x - y));
                        }
                        first_el.extend(rest);
                        *val = first_el;
                    }
                }
            }
            return ();
        }

        fn calc_beta(&mut self) -> () {
            for a in self.entropy_vectors.iter_mut() {
                for (_, val) in a.entropy_values.iter_mut() {
                    if val.len() != 0 {
                        let mut first_el: f64 = val[0];
                        let mut rest: Vec<f64> = Vec::new();
                        rest.push(first_el);
                        for x in val[1..].iter() {
                            first_el = first_el + x;
                            rest.push(first_el);
                        }
                        *val = rest;
                    }
                }
            }
            return ();
        }

        fn calc_q(&mut self) -> () {
            for a in self.entropy_vectors.iter_mut() {
                for (_, val) in a.entropy_values.iter_mut() {
                    if val.len() != 0 {
                        let sum: f64 = val.iter().sum();
                        let mut rest: Vec<f64> = Vec::new();
                        for x in val.iter() {
                            rest.push(x / sum);
                        }
                        *val = rest;
                    }
                }
            }
            return ();
        }

        fn calc_H(&mut self) -> () {
            for a in self.entropy_vectors.iter_mut() {
                for (_, val) in a.entropy_values.iter_mut() {
                    if val.len() != 0 {
                        let len = val.len() as f64;
                        let mut rest: Vec<f64> = Vec::new();
                        for x in val.iter() {
                            rest.push(x * x.log2());
                        }
                        let pre_scale: f64 = -rest.iter().sum::<f64>();
                        if len.log2() != 0.0 {
                            let post_scale: f64 = (pre_scale - 0.0) / (len.log2() - 0.0);
                            *val = vec![post_scale];
                        } else {
                            let post_scale: f64 = 1.0;
                            *val = vec![post_scale];
                        }
                    }
                }
            }
            return ();
        }

        fn chained_operations(&mut self) ->(){
            self.get_position();
            self.calc_alpha();
            self.calc_beta();
            self.calc_q();
            self.calc_H();
            println!("Record {:?} done",self.id)
        }
    }

    #[derive(Clone, Debug)]
    pub struct FeatureExtraction {
        pub origin: String,
        pub position: Vec<u32>,
        pub matrix_dim: u32,
        pub sub_matrix: u32,
        pub full_matrix: Vec<Vec<Array2<f64>>>,
        pub new_matrix: Vec<Vec<Array2<f64>>>,
        pub feature_vector: Vec<Vec<Vec<u32>>>,
    }
    impl FeatureExtraction {
        pub fn new(record: &Record,config: &Generate_Config) -> Result<FeatureExtraction,&'static str>{
            let origin=record.id.clone();
            let matrix_dim=config.matrix_dim.clone();
            let sub_matrix=config.sub_matrix.clone();
            let position=FeatureExtraction::return_position_vector(record.entropy_vectors.len()as u32,config.matrix_dim);
            let new_matrix= match FeatureExtraction::return_empty_new_matrix(&record,config.matrix_dim){Ok(v)=>v,Err(_e)=>return Err("Error generating empty matrix")};
            let feature_vector=match FeatureExtraction::return_empty_feature_vector(&record,config.matrix_dim){Ok(v)=>v,Err(_e)=>return Err("Error generating feature vector")};
            let full_matrix=match FeatureExtraction::return_full_matrix(record,config.matrix_dim){Ok(v)=> v, Err(_e)=> return Err("Error generating full matrix")};
            return Ok(FeatureExtraction{origin,position,matrix_dim,sub_matrix,full_matrix,new_matrix,feature_vector});
        }

        fn return_position_vector(length: u32, mdim : u32) -> Vec<u32>{
            let rest =length%mdim;
            let position_vec=(0..length-rest).step_by(mdim as usize).collect();
            return position_vec;
        }

        fn return_full_matrix(record: &Record, mdim : u32) -> Result<Vec<Vec<Array2<f64>>>,&'static str>{
            let mykmerdict=&record.entropy_vectors[0].kmer_dict.clone();
            let mykmer=&record.entropy_vectors[0].kmers.clone();
            let mut outer=Vec::new();
            for ch in record.entropy_vectors.chunks(mdim as usize){
                if ch.len() as u32 != mdim{continue};
                let mut inner=Vec::new();
                for kmer in mykmer{
                    let mut ordered_kmers=match mykmerdict.get(&kmer.to_string()){Some(v) => v.clone(), None=> return Err("failure retriving from Entropy.kmer_dict")};
                    ordered_kmers.sort();
                    let pre_matrix: Vec<Vec<f64>>=ch.iter().map(|x| ordered_kmers.iter().map(|y| if !x.entropy_values.get(y).unwrap().is_empty(){
                        x.entropy_values.get(y).unwrap()[0]}else{
                         0.0}).collect()).collect();
                    let mat=Array::from_shape_vec((pre_matrix[0].len(),pre_matrix.len()).f(),pre_matrix.into_iter().flatten().collect());
                    inner.push(mat.unwrap());
                }
                outer.push(inner);
            }
            return Ok(outer);
        }

        pub fn return_new_matrix(&mut self)->(){
            for (olfull,olnew) in self.full_matrix.iter().zip(self.new_matrix.iter_mut()){
                for (ilfull,ilnew) in olfull.iter().zip(olnew.iter_mut()){
                    let mean=ilfull.sum();
                    *ilnew=ilfull.clone();
                    ilnew.mapv_inplace(|a| if a <mean{0.0}else{a});
                }
            }
        }

        pub fn return_feature_vector(&mut self){
            for (olnew,olfvec) in self.new_matrix.iter().zip(self.feature_vector.iter_mut()){
                for (ilnew,ilfvec) in olnew.iter().zip(olfvec.iter_mut()){
                    let mut pre_fvec : Vec<f64>= Vec::new();
                    for i in (0..ilnew.nrows()).step_by(self.sub_matrix as usize){
                        for j in (0..ilnew.ncols()).step_by(self.sub_matrix as usize){
                            let var = ilnew.slice(s![i..i+self.sub_matrix as usize,j..j+self.sub_matrix as usize]).sum();
                            pre_fvec.push(var);
                        }
                    }
                    *ilfvec=FeatureExtraction::argsort(pre_fvec);
                } 
            }
        }

        fn argsort(vec:Vec<f64>) -> Vec<u32>{
            let mut sorted_vec:Vec<(usize,f64)>=vec.iter().enumerate().map(|(i, v)| (i,*v)).collect::<Vec<(usize,f64)>>();
            sorted_vec.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
            let out_vec:Vec<u32>=sorted_vec.iter().rev().take(6).map(|a| a.0 as u32).collect();
            return out_vec;
        }

        fn return_empty_new_matrix(record:&Record, mdim : u32)->Result<Vec<Vec<Array2<f64>>>,&'static str>{
            let inner : Vec<Array2<f64>>=vec![Array2::<f64>::zeros((0,0));record.kmers.len()];
            let outer:Vec<Vec<Array2<f64>>>=vec![inner;record.entropy_vectors.len()/mdim as usize];
            return Ok(outer);
        }

        fn return_empty_feature_vector(record:&Record, mdim : u32)->Result<Vec<Vec<Vec<u32>>>,&'static str>{
            let inner : Vec<Vec<u32>>=vec![vec![0];record.kmers.len()];
            let outer:Vec<Vec<Vec<u32>>>=vec![inner;record.entropy_vectors.len()/mdim as usize];
            return Ok(outer);
        }

    }
   
    #[derive(Clone, Debug)]
    pub struct Point{
        abs_position : u32,
        id_red : u32,
        fvec:Vec<u32>,
        couple:ad_couple,
    }
    
    #[derive(Clone, Debug)]
    enum ad_couple{
        db_couple(u32,String),
        query_couple(u32),
    }

    impl Point{
        pub fn from_FeatureExtraction(feature:&FeatureExtraction,config : &Generate_Config)-> Result<Vec<Point>,&'static str> {
            let mut out : Vec<Point> = Vec::new();
            let myposition : Vec<u32>=feature.position.iter().flat_map(|x| std::iter::repeat(*x).take(feature.feature_vector[0].len())).collect();
            for ((id_red,fvec),position) in feature.feature_vector.iter().flatten().enumerate().zip(myposition.iter()){
                out.push(Point::new(id_red as u32,fvec.to_vec(),*position,config,&feature.origin));
            }
            return Ok(out);
        }

        pub fn new(id_red:u32,fvec:Vec<u32>,position:u32,config:&Generate_Config,id:&String)->Point{
            let mycouple = match config.database{true => ad_couple::db_couple(position,id.clone()),false=>ad_couple::query_couple(position)};
            return Point{abs_position:position,id_red:id_red,fvec:fvec,couple:mycouple};
        }


    }

    //functions inside shannon_entropy module

    fn make_adresses(points : Vec<Point> , config : &Generate_Config) ->Vec<(String,ad_couple)>{
            let mut adresses : Vec<(String,ad_couple)> = Vec::new();
            for (zone,anchor) in points[(config.offset as usize)..].windows(5).zip(points.iter()){
                for z in zone.iter(){
                    adresses.push(([anchor.fvec.iter().join(""),z.fvec.iter().join(""),(anchor.abs_position as i32 - z.abs_position as i32).abs().to_string()].join(""),z.couple.clone()))
                }
            } 
            return adresses;
        }

    pub fn entropy_calculation_main(config: Generate_Config) -> Result<(), &'static str> {
        let mut records_in_file = match read_data(&config) {
            Ok(v) => v,
            Err(e) => return Err(e),
        };
        records_in_file
            .par_iter_mut()
            .for_each(|x| x.chained_operations());
        if config.early_exit {
            print_to_file(&records_in_file, &config.output_dir)
        };
        let mut features:Vec<FeatureExtraction>=records_in_file.iter().map(|x| FeatureExtraction::new(x,&config).unwrap()).collect();
        features.par_iter_mut().for_each(|x| x.return_new_matrix());
        features.par_iter_mut().for_each(|x| x.return_feature_vector());
        let points : Vec<Point> = features.iter().map(|x| Point::from_FeatureExtraction(x,&config).unwrap()).collect::<Vec<Vec<Point>>>().into_iter().flatten().collect();
        let adresses : Vec<(String,ad_couple)> = make_adresses(points,&config);
        print_adresses_to_file(adresses,&config.output_dir);
        return Ok(());
    }

    fn print_adresses_to_file(adr: Vec<(String,ad_couple)>,outpath: &str)->(){
        let path = Path::new(outpath).join("dnasr_adresses.txt");
        println!("Printing adresses to :{:?}", path.display());
        let newfile = match File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}", path.display(), why),
            Ok(file) => file,
        };
        let mut writer = BufWriter::new(&newfile);
        for a in adr{
            let c = match a.1{ad_couple::db_couple(x,y) =>(x,y),ad_couple::query_couple(x)=>(x,String::from(""))};
            write!(&mut writer, "{}\t{}\t{}\n",a.0,c.0,c.1);
        }
    }

    fn read_data(config: &Generate_Config) -> Result<Vec<Record>, &'static str> {
        let mut record_vector: Vec<Record> = Vec::new();
        if config.input_file.contains(".fasta") {
            let reader = fasta::Reader::from_file(&config.input_file);
            if let Err(_e) = reader {
                return Err("Error in fasta::Reader");
            };
            let records = reader.unwrap().records().map(|record| record.unwrap());
            for record in records {
                let new_record = match Record::new(config, record.seq(), record.id()) {
                    Ok(w) => w,
                    Err(_e) => return Err("Error generating window struct"),
                };
                record_vector.push(new_record);
            }
        } else if config.input_file.contains(".fastq") {
            let reader = fastq::Reader::from_file(&config.input_file);
            if let Err(_e) = reader {
                return Err("Error in fastq::Reader");
            };
            let records = reader.unwrap().records().map(|record| record.unwrap());
            for record in records {
                let new_record = match Record::new(config, record.seq(), record.id()) {
                    Ok(w) => w,
                    Err(_e) => return Err("Error generating window struct"),
                };
                record_vector.push(new_record);
            }
        }
        return Ok(record_vector);
    }

    fn print_to_file(records: &Vec<Record>, outpath: &str) -> () {
        let path = Path::new(outpath).join("dnasr_early_exit.txt");
        println!("Printing early exit content to :{:?}", path.display());
        let newfile = match File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}", path.display(), why),
            Ok(file) => file,
        };
        let mut writer = BufWriter::new(&newfile);
        for rec in records.iter() {
            for (win, ent) in
                (1..rec.window_starts.iter().count() as u32).zip(rec.entropy_vectors.iter())
            {
                for (key, val) in ent.entropy_values.iter() {
                    let out = match val.len() {
                        0 => 0.0,
                        _ => val[0],
                    };
                    write!(&mut writer, "{}\t{}\t{}\t{}\n", rec.id, win, key, out)
                        .expect("Not written");
                }
            }
        }
        println!("Exiting....");
        process::exit(0);
    }
}

pub mod search {
    use super::*;
    pub fn search_main(config: Search_Config) -> Result<(), &'static str> {
        return Err("not implemented yet");
    }
}
