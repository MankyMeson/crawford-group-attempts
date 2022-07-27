use std::io;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Error;
use std::io::ErrorKind;
use std::env;
use std::fs;

fn main() -> io::Result<()> {

    let mut args = env::args();
    args.next();

    for arg in args {
        let lines = file_to_vec(arg)?;

        // Read first line (No. atoms in system)
        let natoms: i32 = (&lines[0]).parse().unwrap();

        if (lines.len()-1) as i32 != natoms {
            return Err(
                Error::new(
                    ErrorKind::InvalidData,
                    "Number of atoms != number of data points"
                )
            );
        }

        // Read all other lines to array of atom data
        // negative Z values are used as errors
        let mut ions: Vec<Ion> = vec![
            Ion {z_val:-1_i32,x:0.,y:0.,z:0.};
            lines.len()-1
        ];

        for i in 1..lines.len() {
            let mut ion_data = (&lines[i]).split_whitespace();

            ions[i-1] = Ion {
                z_val: ion_data.next().unwrap().parse::<i32>().unwrap(),
                x: ion_data.next().unwrap().parse::<f64>().unwrap(),
                y: ion_data.next().unwrap().parse::<f64>().unwrap(),
                z: ion_data.next().unwrap().parse::<f64>().unwrap(),
            }

        }

        let bond_lengths = all_bond_lengths(&ions).unwrap();
        let bond_angles = bond_angles(&ions).unwrap();

        //testing output
        println!("number of atoms:\n    {:?}", natoms);
        println!("ion data:\n   {:?}", ions);
        println!("all bond lengths:\n    {:?}", bond_lengths);
        println!("all bond angles:\n    {:?}", bond_angles);
    }

    Ok(())
}


fn file_to_vec(filename: String) -> io::Result<Vec<String>> {
    let file_in = fs::File::open(filename)?;
    let file_reader = BufReader::new(file_in);

    Ok(file_reader.lines().filter_map(io::Result::ok).collect())
}


#[derive(Debug,Clone,Copy,PartialEq)]
pub struct Ion {
    pub z_val: i32,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl Ion {

    pub fn bond_length(&self, other: &Ion) -> f64 {
        ( (self.x-other.x).powi(2)
        + (self.y-other.y).powi(2)
        + (self.z-other.z).powi(2)
        ).sqrt()
    }

    pub fn bond_angle(ioni: &Ion, ionj: &Ion, ionk: &Ion) -> f64 {
        let norm_const = ioni.bond_length(&ionj);
        let x_ji = ioni.x - ionj.x;
        let y_ji = ioni.y - ionj.y;
        let z_ji = ioni.z - ionj.z;
        let x_jk = ionk.x - ionj.x;
        let y_jk = ionk.y - ionj.y;
        let z_jk = ionk.z - ionj.z;
        ((x_ji*x_jk + y_ji*y_jk + z_ji*z_jk)
        /(norm_const*norm_const)).acos()
    }
}

pub fn all_bond_lengths(mol: &Vec<Ion>) -> io::Result<Vec<Vec<f64>>> {
    if mol.len() <= 1 {
        Err(Error::new(ErrorKind::InvalidData, "No bonds in molecule"))

    } else {
        let mut lengths = vec![vec![0.;mol.len()];mol.len()];

        for i in 0..mol.len() {
            for j in 0..i {
                if i != j {
                    let l = mol[i].bond_length(&mol[j]);
                    lengths[i][j] = l;
                    lengths[j][i] = l;
                }
            }
        }

        Ok(lengths)
    }
}

pub fn bond_angles(mol: &Vec<Ion>) -> io::Result<Vec<(usize,usize,usize,f64)>> {
    let len = mol.len();
    if len <= 2 {
        Err(Error::new(ErrorKind::InvalidData, "too few ions for bond angles"))
    } else {
        // length of angles is nth trigonal pyramidal number
        // where n is len-1
        let n_uniques = (3*(len-2).pow(2) + (len-2).pow(3) + 2*(len-2))/6;
        let mut angles = vec![(0,0,0,0.);n_uniques];
        let mut angle_index = 0;
        for i in 0..mol.len() {
            for j in 0..i {
                for k in 0..j {
                    if !(i==j||j==k||i==k) {
                        angles[angle_index] = (k,j,i,Ion::bond_angle(&mol[i],&mol[j],&mol[k]));
                        angle_index += 1;
                    }
                }
            }
        }
        Ok(angles)
    }
}
