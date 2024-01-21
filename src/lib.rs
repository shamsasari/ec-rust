#![deny(unused_must_use)]

use std::ops::Add;
use num_bigint::{BigUint, ToBigUint};

use crate::Point::{Coordinate, Identity};

#[derive(PartialEq, Debug, Clone)]
enum Point {
    Coordinate(BigUint, BigUint),
    Identity
}

impl Point {
    fn is_on_curve(&self, curve: &EllipticCurve) -> bool {
        match self {
            Coordinate(x, y) => {
                // y^2 = x^3 + ax + b equation of the elliptic curve
                let x3 = x.mod_pow(3, &curve.p);
                let ax = &curve.a * x;
                y.square(&curve.p) == (x3 + ax + &curve.b) % &curve.p
            },
            Identity => true
        }
    }
}

struct EllipticCurve {
    // y^2 = x^2 + ax + b
    a: BigUint,
    b: BigUint,
    p: BigUint,
}

impl EllipticCurve {
    fn add(&self, point1: &Point, point2: &Point) -> Point {
        assert!(point1.is_on_curve(self), "point1 is not on the curve");
        assert!(point2.is_on_curve(self), "point2 is not on the curve");

        match (point1, point2) {
            (Identity, _) => point2.clone(),
            (_, Identity) => point1.clone(),
            (Coordinate(x1, y1), Coordinate(x2, y2)) => {
                let s = if point1 != point2 {
                    div(&sub(&y2, &y1, &self.p), &sub(&x2, &x1, &self.p), &self.p)
                } else {
                    (3u32 * x1.square(&self.p) - &self.a) / 2u32 * y1
                };
                // x3 = s^2 - x1 - x2 mod p
                let x3 = sub(&sub(&s.square(&self.p), &x1, &self.p), &x2, &self.p);
                // y3 = s(x1 - x3) - y1 mod p
                let y3 = sub(&mult(&s, &sub(&x1, &x3, &self.p), &self.p), &y1, &self.p);
                Coordinate(x3, y3)
            }
        }
    }

    fn double(point: &Point) -> Point {
        todo!()
    }

    fn scalar_mult(point: &Point, int: &BigUint) -> Point {
        todo!()
    }
}

struct FiniteField {

}

impl FiniteField {

}

fn add(int1: &BigUint, int2: &BigUint, p: &BigUint) -> BigUint {
    (int1 + int2) % p
}

fn mult(int1: &BigUint, int2: &BigUint, p: &BigUint) -> BigUint {
    (int1 * int2) % p
}

fn sub(int1: &BigUint, int2: &BigUint, p: &BigUint) -> BigUint {
    (int1 + int2.additive_inverse(&p)) % p
}

fn div(int1: &BigUint, int2: &BigUint, p: &BigUint) -> BigUint {
    mult(int1, &int2.multiplicative_inverse(&p), &p)
}

#[derive(Debug)]
struct ModInt {
    value: BigUint,
    modulus: BigUint
}

impl ModInt {
    fn of(value: i32, modulus: i32) -> Self {
        ModInt { value: big_uint(value), modulus: big_uint(modulus) }
    }

    // fn additive_inverse(&self, modulus: &BigUint) -> BigUint {
    //     assert!(self < modulus, "modulus must be larger");
    //     modulus - self
    // }
}

impl Add<ModInt> for ModInt {
    type Output = Self;

    fn add(self, rhs: ModInt) -> Self::Output {
        assert_eq!(self.modulus, rhs.modulus);
        let modulus = self.modulus.clone();
        ModInt { value: (self.value + rhs.value) % &modulus, modulus }
    }
}

trait BigUintExt {
    fn additive_inverse(&self, modulus: &BigUint) -> BigUint;
    fn multiplicative_inverse(&self, p: &BigUint) -> BigUint;
    fn square(&self, modulus: &BigUint) -> BigUint;
    fn mod_pow(&self, exponent: i32, modulus: &BigUint) -> BigUint;
}

impl BigUintExt for BigUint {
    fn additive_inverse(&self, modulus: &BigUint) -> BigUint {
        assert!(self < modulus, "modulus must be larger");
        modulus - self
    }

    fn multiplicative_inverse(&self, p: &BigUint) -> BigUint {
        // TODO check p is prime
        // int^(-1) mod p = int^(p-2) mod p (fermat's little theorem)
        self.modpow(&(p - 2u32), p)
    }

    fn square(&self, modulus: &BigUint) -> BigUint {
        self.modpow(&big_uint(2), modulus)
    }
    fn mod_pow(&self, exponent: i32, modulus: &BigUint) -> BigUint {
        self.modpow(&big_uint(exponent), modulus)
    }
}

fn big_uint(value: i32) -> BigUint {
    BigUint::from(value as u32)
}

#[cfg(test)]
mod tests {
    use num_traits::identities::Zero;

    use super::*;

    #[test]
    fn test_add() {
        let a = big_uint(5);
        let b = big_uint(10);
        let p = big_uint(11);
        assert_eq!(add(&a, &b, &p), big_uint(4));
    }

    #[test]
    fn test_mult() {
        let a = big_uint(5);
        let b = big_uint(10);
        let p = big_uint(11);
        assert_eq!(mult(&a, &b, &p), big_uint(6));
    }

    #[test]
    fn additive_inverse() {
        let a = big_uint(9);
        let p = big_uint(11);
        assert_eq!(a.additive_inverse(&p), big_uint(2));
    }

    #[test]
    fn additive_inverse_to_get_identity() {
        let a = big_uint(9);
        let p = big_uint(11);
        assert_eq!(add(&a, &a.additive_inverse(&p), &p), BigUint::zero());
    }

    #[test]
    fn subtract_1() {
        let a = big_uint(4);
        let b = big_uint(9);
        let p = big_uint(11);
        assert_eq!(sub(&a, &b, &p), big_uint(6));
    }

    #[test]
    fn subtract_2() {
        let a = big_uint(15);
        let b = big_uint(2);
        let p = big_uint(11);
        assert_eq!(sub(&a, &b, &p), big_uint(2));
    }

    #[test]
    fn multiplicative_inverse_to_get_identity() {
        let a = big_uint(4);
        let p = big_uint(11);
        assert_eq!(mult(&a, &a.multiplicative_inverse(&p), &p), big_uint(1));
    }

    #[test]
    fn ec_point_addition() {
        // y^2 = x^3 + 2x + 2 mod 17
        let ec = EllipticCurve {
            a: big_uint(2),
            b: big_uint(2),
            p: big_uint(17)
        };

        // (6, 3) + (5, 1) = (10, 6)

        let point1 = Coordinate(big_uint(6), big_uint(3));
        let point2 = Coordinate(big_uint(5), big_uint(1));

        assert_eq!(ec.add(&point1, &point2), Coordinate(big_uint(10), big_uint(6)));
    }
}

