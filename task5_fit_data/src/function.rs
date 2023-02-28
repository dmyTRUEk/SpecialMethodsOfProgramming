//! Function struct.

use rand::{Rng, thread_rng};

use crate::{
    extensions::ExtGenFromArray,
    float_type::float,
    param::{PARAMETER_NAMES, ParamName},
    params::Params,
};


#[derive(Debug, Clone, PartialEq)]
pub enum Function {
    X,
    Const { value: float }, Zero, One,
    Param { name: ParamName },

    Neg { value: Box<Self> },

    Exp { value: Box<Self> },
    Ln  { value: Box<Self> }, // natural
    Sqrt{ value: Box<Self> }, // square root
    Sq  { value: Box<Self> }, // squared

    Sin { value: Box<Self> },
    Cos { value: Box<Self> },
    Tan { value: Box<Self> },

    // TODO(feat): add more functions.

    // ArcSin { value: Box<Self> },
    // ArcCos { value: Box<Self> },
    // ArcTan { value: Box<Self> },

    // Sinh { value: Box<Self> },
    // Cosh { value: Box<Self> },
    // Tanh { value: Box<Self> },

    Add { lhs: Box<Self>, rhs: Box<Self> },
    Sub { lhs: Box<Self>, rhs: Box<Self> },
    Mul { lhs: Box<Self>, rhs: Box<Self> },
    Div { lhs: Box<Self>, rhs: Box<Self> },
    Pow { lhs: Box<Self>, rhs: Box<Self> },

    // special cases:

    #[allow(unused)]
    /// a + bx + cx^2 + dx^3 + …
    Polynomial { degree: usize },

    #[allow(unused)]
    /// a + bx + cx^2/2! + dx^3/3! + …
    BtrPolynomial { degree: usize },

    // TODO
    // FourierConstSinCosSeries { degree: usize },
}

impl Function {
    pub fn from_str(string: &str) -> Result<Self, String> {
        const DEBUG: bool = false;
        const LB: char = '(';
        const RB: char = ')';
        const PLUS: char = '+';
        const MINUS: char = '-';
        const MULTIPLY: char = '*';
        const DIVIDE: char = '/';
        const POWER: char = '^';
        if DEBUG { println!("{}", string) }
        let string: String = string.replace(' ', "");
        let len = string.len();

        if string.chars().next().unwrap() == LB && string.chars().last().unwrap() == RB {
            // check if first LB corresponds to last RB
            let mut is_first_lb_corresponds_to_last_rb: bool = true;
            let mut brackets_level: i32 = 0;
            for c in string[..len-1].chars() {
                match c {
                    LB => brackets_level += 1,
                    RB => brackets_level -= 1,
                    _ => {}
                }
                if brackets_level == 0 {
                    is_first_lb_corresponds_to_last_rb = false;
                    break;
                }
            }
            if is_first_lb_corresponds_to_last_rb {
                return Self::from_str(&string[1..len-1]);
            }
        }

        if string.contains([LB, RB, PLUS, MINUS, MULTIPLY, DIVIDE, POWER]) {
            if string.matches(LB).count() != string.matches(RB).count() { return Err("Bad brackets sequence: opening and closing brackets count doesn't match.".to_string()) }
            let chars = string.chars().collect::<Vec<_>>(); // this is needed for reverse iterator.
            // check brackets level
            let mut brackets_level: i32 = 0;
            for &c in chars.iter() {
                match c {
                    LB => brackets_level += 1,
                    RB => brackets_level -= 1,
                    _ => {}
                }
                if brackets_level < 0 { return Err("Bad brackets sequence: closing bracket found before opening one.".to_string()) }
            }
            // "split" by PLUS
            let mut brackets_level: i32 = 0;
            for (i, &c) in chars.iter().enumerate() {
                match c {
                    LB => brackets_level += 1,
                    RB => brackets_level -= 1,
                    PLUS if brackets_level == 0 => {
                        return Ok(Self::Add {
                            lhs: box Self::from_str(&string[..i]).unwrap(),
                            rhs: box Self::from_str(&string[i+1..]).unwrap()
                        });
                    }
                    _ => {}
                }
            }
            // "split" by MINUS
            let mut brackets_level: i32 = 0;
            for (i, &c) in chars.iter().enumerate().rev() {
                match c {
                    LB => brackets_level += 1,
                    RB => brackets_level -= 1,
                    MINUS if brackets_level == 0 && i != 0 => {
                        return Ok(Self::Sub {
                            lhs: box Self::from_str(&string[..i]).unwrap(),
                            rhs: box Self::from_str(&string[i+1..]).unwrap()
                        });
                    }
                    _ => {}
                }
            }

            // NEG must be before POWER (for `-x^2` to work) but after MINUS (for `-a+b` to work).
            if string.starts_with(MINUS) { return Ok(Self::Neg { value: box Self::from_str(&string[1..])? }) }

            // "split" by MULTIPLY
            let mut brackets_level: i32 = 0;
            for (i, &c) in chars.iter().enumerate() {
                match c {
                    LB => brackets_level += 1,
                    RB => brackets_level -= 1,
                    MULTIPLY if brackets_level == 0 => {
                        return Ok(Self::Mul {
                            lhs: box Self::from_str(&string[..i]).unwrap(),
                            rhs: box Self::from_str(&string[i+1..]).unwrap()
                        });
                    }
                    _ => {}
                }
            }
            // "split" by DIVIDE
            let mut brackets_level: i32 = 0;
            for (i, &c) in chars.iter().enumerate() {
                match c {
                    LB => brackets_level += 1,
                    RB => brackets_level -= 1,
                    DIVIDE if brackets_level == 0 => {
                        return Ok(Self::Div {
                            lhs: box Self::from_str(&string[..i]).unwrap(),
                            rhs: box Self::from_str(&string[i+1..]).unwrap()
                        });
                    }
                    _ => {}
                }
            }
            // "split" by POW
            let mut brackets_level: i32 = 0;
            for (i, c) in string.chars().enumerate() {
                match c {
                    LB => brackets_level += 1,
                    RB => brackets_level -= 1,
                    POWER if brackets_level == 0 => {
                        let lhs = box Self::from_str(&string[..i]).unwrap();
                        let rhs = box Self::from_str(&string[i+1..]).unwrap();
                        return Ok(
                            if *rhs == (Function::Const { value: 2. }) {
                                Self::Sq { value: lhs }
                            } else {
                                Self::Pow { lhs, rhs }
                            }
                        );
                    }
                    _ => {}
                }
            }
        }
        else {
            if string == "x" { return Ok(Function::X) }
            let first_char: char = string.chars().next().unwrap();
            if len == 1 && PARAMETER_NAMES.contains(&first_char) {
                return Ok(Function::Param { name: first_char });
            }
            if string == "0" { return Ok(Function::Zero) }
            if string == "1" { return Ok(Function::One) }
            match string.parse::<float>() {
                Ok(number) => return Ok(Function::Const { value: number }),
                Err(_e) => {}
            }
        }

        if string.starts_with("exp(") && string.ends_with(RB) { return Ok(Self::Exp { value: box Self::from_str(&string[4..len-1])? }) }
        if string.starts_with("ln(")  && string.ends_with(RB) { return Ok(Self::Ln  { value: box Self::from_str(&string[3..len-1])? }) }
        if string.starts_with("sqrt(")&& string.ends_with(RB) { return Ok(Self::Sqrt{ value: box Self::from_str(&string[5..len-1])? }) }
        if string.ends_with("^2")                             { return Ok(Self::Sq  { value: box Self::from_str(&string[0..len-2])? }) }
        if string.starts_with("sin(") && string.ends_with(RB) { return Ok(Self::Sin { value: box Self::from_str(&string[4..len-1])? }) }
        if string.starts_with("cos(") && string.ends_with(RB) { return Ok(Self::Cos { value: box Self::from_str(&string[4..len-1])? }) }
        if string.starts_with("tan(") && string.ends_with(RB) { return Ok(Self::Tan { value: box Self::from_str(&string[4..len-1])? }) }

        Err("Unable to parse.".to_string())
    }


    pub fn gen(complexity: u32) -> Self {
        // TODO?: add arg `non_param` which disallows it to be param (for unary and one side of binary?).
        let mut rng = thread_rng();
        match complexity {
            0 => {
                match rng.gen_range(0 ..= 1) {
                    0 => Self::X,
                    1 => Self::Param { name: rng.gen_from_array(PARAMETER_NAMES) },
                    // 2 => Self::Const {
                    //     value: if rng.gen_range(0 .. 100) < 50 { 2. } else { rng.gen_range(-3. ..= 5.) }
                    // },
                    _ => unreachable!()
                }
            }
            _ => {
                let complexity = complexity - 1;
                let partition = rng.gen_range(0 ..= complexity);
                match rng.gen_range(0 ..= 12) {
                    0 => Self::Neg { value: box Self::gen(complexity) },

                    1 => Self::Exp { value: box Self::gen(complexity) },
                    2 => Self::Ln  { value: box Self::gen(complexity) },
                    3 => Self::Sqrt{ value: box Self::gen(complexity) },
                    4 => Self::Sq  { value: box Self::gen(complexity) },

                    5 => Self::Sin { value: box Self::gen(complexity) },
                    6 => Self::Cos { value: box Self::gen(complexity) },
                    7 => Self::Tan { value: box Self::gen(complexity) },

                    8 => Self::Add { lhs: box Self::gen(partition), rhs: box Self::gen(complexity-partition) },
                    9 => Self::Sub { lhs: box Self::gen(partition), rhs: box Self::gen(complexity-partition) },
                    10 => Self::Mul { lhs: box Self::gen(partition), rhs: box Self::gen(complexity-partition) },
                    11 => Self::Div { lhs: box Self::gen(partition), rhs: box Self::gen(complexity-partition) },
                    12 => Self::Pow { lhs: box Self::gen(partition), rhs: box Self::gen(complexity-partition) },
                    _ => unreachable!()
                }
            }
        }
    }


    // #[inline(always)]
    // TODO(optimize)?: `pub const fn`
    pub fn eval(&self, x: float, params: &Params) -> float {
        const DEBUG: bool = false;
        let result = match self {
            Self::X => x,
            Self::Const { value } => *value,
            Self::Zero => 0.,
            Self::One  => 1.,
            Self::Param { name } => params.get_by_name_unchecked(*name),

            Self::Neg { value } => -value.eval(x, params),

            Self::Exp { value } => value.eval(x, params).exp(),
            Self::Ln  { value } => value.eval(x, params).ln(),
            Self::Sqrt{ value } => value.eval(x, params).sqrt(),
            Self::Sq  { value } => value.eval(x, params).powi(2),

            Self::Sin { value } => value.eval(x, params).sin(),
            Self::Cos { value } => value.eval(x, params).cos(),
            Self::Tan { value } => value.eval(x, params).tan(),

            Self::Add { lhs, rhs } => lhs.eval(x, params) + rhs.eval(x, params),
            Self::Sub { lhs, rhs } => lhs.eval(x, params) - rhs.eval(x, params),
            Self::Mul { lhs, rhs } => lhs.eval(x, params) * rhs.eval(x, params),
            Self::Div { lhs, rhs } => lhs.eval(x, params) / rhs.eval(x, params),
            Self::Pow { lhs, rhs } => lhs.eval(x, params).powf(rhs.eval(x, params)),

            Self::Polynomial { degree } => {
                let mut r: float = 0.;
                for i in 0..=*degree {
                    r += params.get_by_name_unchecked(PARAMETER_NAMES[i]) * x.powi(i as i32);
                }
                r
            }
            Self::BtrPolynomial { degree } => {
                let mut r: float = 0.;
                let mut denominator: float = 1.;
                for i in 0..=*degree {
                    if i >= 2 {
                        denominator *= i as float;
                    }
                    r += params.get_by_name_unchecked(PARAMETER_NAMES[i]) * x.powi(i as i32) / denominator;
                }
                r
            }
        };
        if DEBUG {
            println!("Function::eval: f(x) = {}\tx={}\t{:?}", self.to_string(), x, params);
            println!("Function::eval: result = {}", result);
        }
        result
    }


    #[allow(unused)]
    // #[inline(always)]
    pub fn eval2(&self, x: float, params: &Params) -> float {
        const STACK_INIT_SIZE: usize = 20;
        type Operand = float;
        // #[repr(u8)] // TODO: try enable for better perfomance
        #[derive(Debug)]
        enum Operator {
            LeftBracket, RightBracket,
            Plus, Minus, Multiply, Divide, Power,
            Neg, Exp, Ln, Sqrt, Sq, Sin, Cos, Tan,
            Polynomial, BtrPolynomial,
        }
        #[derive(Debug)]
        enum Token {
            Operand(Operand),
            Operator(Operator)
        }
        fn gen_token_list(f: &Function, x: float, p: &Params) -> Vec<Token> {
            let self_: Token = match f {
                Function::X => Token::Operand(x),
                Function::Const { value } => Token::Operand(*value),
                Function::Zero => Token::Operand(0.),
                Function::One  => Token::Operand(1.),
                Function::Param { name } => Token::Operand(p.get_by_name_unchecked(*name)),

                Function::Neg { .. } => Token::Operator(Operator::Neg),

                Function::Exp { .. } => Token::Operator(Operator::Exp),
                Function::Ln  { .. } => Token::Operator(Operator::Ln),
                Function::Sqrt{ .. } => Token::Operator(Operator::Sqrt),
                Function::Sq  { .. } => Token::Operator(Operator::Sq),

                Function::Sin { .. } => Token::Operator(Operator::Sin),
                Function::Cos { .. } => Token::Operator(Operator::Cos),
                Function::Tan { .. } => Token::Operator(Operator::Tan),

                Function::Add { .. } => Token::Operator(Operator::Plus),
                Function::Sub { .. } => Token::Operator(Operator::Minus),
                Function::Mul { .. } => Token::Operator(Operator::Multiply),
                Function::Div { .. } => Token::Operator(Operator::Divide),
                Function::Pow { .. } => Token::Operator(Operator::Power),
                _ => todo!()
            };
            let mut children: Vec<Token> = vec![];
            // TODO?: `Vec<Token>` -> Option<Vec<Token>>
            match f {
                Function::X
                | Function::Const { .. }
                | Function::Zero
                | Function::One
                | Function::Param { .. }
                => {},

                Function::Neg { value }
                | Function::Exp { value }
                | Function::Ln  { value }
                | Function::Sqrt{ value }
                | Function::Sq  { value }
                | Function::Sin { value }
                | Function::Cos { value }
                | Function::Tan { value }
                => {
                    children.push(Token::Operator(Operator::LeftBracket));
                    children.extend(gen_token_list(value, x, p));
                    children.push(Token::Operator(Operator::RightBracket));
                }

                Function::Add { lhs, rhs }
                | Function::Sub { lhs, rhs }
                | Function::Mul { lhs, rhs }
                | Function::Div { lhs, rhs }
                | Function::Pow { lhs, rhs }
                => {
                    children.push(Token::Operator(Operator::LeftBracket));
                    children.extend(gen_token_list(lhs, x, p));
                    children.extend(gen_token_list(rhs, x, p));
                    children.push(Token::Operator(Operator::RightBracket));
                }
                _ => todo!()
            };
            enum TraversalType {
                PreOrder, InOrder, PostOrder,
                PreOrderRev, InOrderRev, PostOrderRev,
            }
            const TRAVERSAL_TYPE: TraversalType = TraversalType::PreOrder;
            let mut tokens: Vec<Token> = match TRAVERSAL_TYPE {
                // TODO(optimize)?: replace `.chain()` with `Vec::append`, `Vec::extend` or smt like that.
                TraversalType::PreOrder  | TraversalType::PreOrderRev  => vec![self_].into_iter().chain(children.into_iter()).collect(),
                TraversalType::PostOrder | TraversalType::PostOrderRev => children.into_iter().chain(vec![self_].into_iter()).collect(),
                _ => todo!()
            };
            match TRAVERSAL_TYPE {
                TraversalType::PreOrderRev
                | TraversalType::InOrderRev
                | TraversalType::PostOrderRev
                => tokens.reverse(),
                _ => {}
            }
            tokens
        }
        println!("{}", self.to_string());
        let mut tokens = gen_token_list(self, x, params);
        println!("{:#?}", tokens);

        let mut operand_stack = Vec::<Operand>::with_capacity(STACK_INIT_SIZE);
        let mut operator_stack = Vec::<Operator>::with_capacity(STACK_INIT_SIZE);

        for token in tokens {
            match token {
                Token::Operand(operand) => operand_stack.push(operand),
                Token::Operator(operator) => {
                    match operator {
                        // Operator::RightBracket => 
                        _ => operator_stack.push(operator)
                    }
                }
            }
        }

        todo!()
    }


    pub fn get_params_names(&self) -> Vec<ParamName> {
        match self {
            Self::X
            | Self::Const { .. }
            | Self::Zero
            | Self::One
            => vec![],

            Self::Param { name } => vec![*name],

            Self::Neg { value }
            | Self::Exp { value }
            | Self::Ln  { value }
            | Self::Sqrt{ value }
            | Self::Sq  { value }
            | Self::Sin { value }
            | Self::Cos { value }
            | Self::Tan { value }
            => value.get_params_names(),

            Self::Add { lhs, rhs }
            | Self::Sub { lhs, rhs }
            | Self::Mul { lhs, rhs }
            | Self::Div { lhs, rhs }
            | Self::Pow { lhs, rhs }
            => [lhs.get_params_names(), rhs.get_params_names()].concat(),

            Self::Polynomial { degree }
            | Self::BtrPolynomial { degree }
            => PARAMETER_NAMES[..=*degree].to_vec(),
        }
    }


    pub fn simplify(self) -> Self {
        const DEBUG: bool = false;
        let mut new_f = self;
        if DEBUG { println!("simplify::begin: {:?}", new_f) }
        new_f = match new_f {
            // not recursive:
            s @ Self::X
            | s @ Self::Const { .. }
            | s @ Self::Zero
            | s @ Self::One
            | s @ Self::Param { .. }
            | s @ Self::Polynomial { .. }
            | s @ Self::BtrPolynomial { .. }
            => s,

            // recursive:
            Self::Neg { value } => Self::Neg { value: box value.simplify() },

            Self::Exp { value } => Self::Exp { value: box value.simplify() },
            Self::Ln  { value } => Self::Ln  { value: box value.simplify() },
            Self::Sqrt{ value } => Self::Sqrt{ value: box value.simplify() },
            Self::Sq  { value } => Self::Sq  { value: box value.simplify() },

            Self::Sin { value } => Self::Sin { value: box value.simplify() },
            Self::Cos { value } => Self::Cos { value: box value.simplify() },
            Self::Tan { value } => Self::Tan { value: box value.simplify() },

            Self::Add { lhs, rhs } => Self::Add { lhs: box lhs.simplify(), rhs: box rhs.simplify() },
            Self::Sub { lhs, rhs } => Self::Sub { lhs: box lhs.simplify(), rhs: box rhs.simplify() },
            Self::Mul { lhs, rhs } => Self::Mul { lhs: box lhs.simplify(), rhs: box rhs.simplify() },
            Self::Div { lhs, rhs } => Self::Div { lhs: box lhs.simplify(), rhs: box rhs.simplify() },
            Self::Pow { lhs, rhs } => Self::Pow { lhs: box lhs.simplify(), rhs: box rhs.simplify() },
        };
        if DEBUG { println!("simplify::middle: {:?}", new_f) }
        new_f = match new_f {
            // simplify params:
            Self::Neg { value: box p @ Self::Param { .. } }   // -a -> a
            | Self::Exp { value: box p @ Self::Param { .. } } // e^a -> a
            | Self::Ln  { value: box p @ Self::Param { .. } } // ln(a) -> a
            | Self::Sqrt{ value: box p @ Self::Param { .. } } // sqrt(a) -> a
            | Self::Sq  { value: box p @ Self::Param { .. } } // sq(a) -> a
            | Self::Sin { value: box p @ Self::Param { .. } } // sin(a) -> a
            | Self::Cos { value: box p @ Self::Param { .. } } // cos(a) -> a
            | Self::Tan { value: box p @ Self::Param { .. } } // tan(a) -> a
            | Self::Add { lhs: box p @ Self::Param { .. }, rhs: box Self::Param { .. } | box Self::Const { .. } | box Self::Zero | box Self::One } // a + b -> a
            | Self::Sub { lhs: box p @ Self::Param { .. }, rhs: box Self::Param { .. } | box Self::Const { .. } | box Self::Zero | box Self::One } // a - b -> a
            | Self::Mul { lhs: box p @ Self::Param { .. }, rhs: box Self::Param { .. } | box Self::Const { .. } | box Self::Zero | box Self::One } // a * b -> a
            | Self::Div { lhs: box p @ Self::Param { .. }, rhs: box Self::Param { .. } | box Self::Const { .. } | box Self::Zero | box Self::One } // a / b -> a
            | Self::Pow { lhs: box p @ Self::Param { .. }, rhs: box Self::Param { .. } | box Self::Const { .. } | box Self::Zero | box Self::One } // a ^ b -> a
            => p,

            Self::Div { lhs, rhs: p @ box Self::Param { .. } } // expr / a -> expr * a
            => Self::Mul { lhs, rhs: p },

            // simplify consts:
            expr @ Self::Neg { value: box Self::Const { .. } | box Self::One | box Self::Zero }   // eval -const
            | expr @ Self::Exp { value: box Self::Const { .. } | box Self::One | box Self::Zero } // eval e^const
            | expr @ Self::Ln  { value: box Self::Const { .. } | box Self::One | box Self::Zero } // eval ln(const)
            | expr @ Self::Sqrt{ value: box Self::Const { .. } | box Self::One | box Self::Zero } // eval e^const
            | expr @ Self::Sq  { value: box Self::Const { .. } | box Self::One | box Self::Zero } // eval e^const
            | expr @ Self::Sin { value: box Self::Const { .. } | box Self::One | box Self::Zero } // eval sin(const)
            | expr @ Self::Cos { value: box Self::Const { .. } | box Self::One | box Self::Zero } // eval cos(const)
            | expr @ Self::Tan { value: box Self::Const { .. } | box Self::One | box Self::Zero } // eval tan(const)
            | expr @ Self::Add { lhs: box Self::Const { .. } | box Self::One | box Self::Zero, rhs: box Self::Const { .. } | box Self::One | box Self::Zero } // eval const + const
            | expr @ Self::Sub { lhs: box Self::Const { .. } | box Self::One | box Self::Zero, rhs: box Self::Const { .. } | box Self::One | box Self::Zero } // eval const - const
            | expr @ Self::Mul { lhs: box Self::Const { .. } | box Self::One | box Self::Zero, rhs: box Self::Const { .. } | box Self::One | box Self::Zero } // eval const * const
            | expr @ Self::Div { lhs: box Self::Const { .. } | box Self::One | box Self::Zero, rhs: box Self::Const { .. } | box Self::One | box Self::Zero } // eval const / const
            | expr @ Self::Pow { lhs: box Self::Const { .. } | box Self::One | box Self::Zero, rhs: box Self::Const { .. } | box Self::One | box Self::Zero } // eval const ^ const
            => Self::Const { value: expr.eval(0., &Params::new()) },

            // simplifies to `0`:
            Self::Sub { lhs: box Self::X, rhs: box Self::X } // x - x == 0
            | Self::Add { lhs: box Self::X, rhs: box Self::Neg { value: box Self::X } } // x + -x == 0
            | Self::Add { lhs: box Self::Neg { value: box Self::X }, rhs: box Self::X } // -x + x == 0
            | Self::Mul { lhs: box Self::Zero, .. } // 0 * expr == 0
            | Self::Mul { rhs: box Self::Zero, .. } // expr * 0 == 0
            | Self::Pow { lhs: box Self::Zero, .. } // 0 ^ expr
            => Self::Zero,

            // simplifies to `1`:
            Self::Div { lhs: box Self::X, rhs: box Self::X } // x / x
            | Self::Div { lhs: box Self::X, rhs: box Self::Neg { value: box Self::X } } // x / -x
            | Self::Div { lhs: box Self::Neg { value: box Self::X }, rhs: box Self::X } // -x / x
            | Self::Pow { lhs: box Self::One, .. }  // 1 ^ expr
            | Self::Pow { rhs: box Self::Zero, .. } // expr ^ 0
            => Self::One,

            // simplify inverse functions: `f(g(expr)) ~= expr`:
            Self::Neg { value: box Self::Neg { value: expr } }   // --expr        == expr
            | Self::Ln  { value: box Self::Exp { value: expr } } // ln(exp(expr)) == expr
            | Self::Exp { value: box Self::Ln  { value: expr } } // exp(ln(expr)) ~= expr
            | Self::Sqrt{ value: box Self::Sq  { value: expr } } // sqrt(sq(expr)) ~= expr
            | Self::Sq  { value: box Self::Sqrt{ value: expr } } // sq(sqrt(expr)) ~= expr
            => *expr,

            // simplifies to `expr`:
            Self::Add { lhs: expr, rhs: box Self::Zero }   // expr + 0
            | Self::Add { lhs: box Self::Zero, rhs: expr } // 0 + expr
            | Self::Sub { lhs: expr, rhs: box Self::Zero } // expr - 0
            | Self::Mul { lhs: expr, rhs: box Self::One }  // expr * 1
            | Self::Mul { lhs: box Self::One, rhs: expr }  // 1 * expr
            | Self::Div { lhs: expr, rhs: box Self::One }  // expr / 1
            | Self::Pow { lhs: expr, rhs: box Self::One }  // expr ^ 1
            => *expr,

            Self::Sub { lhs: box Self::Zero, rhs: expr } // 0 - expr -> -expr
            => Self::Neg { value: expr },

            Self::Sq { value: box Self::Neg { value: expr } } // (-expr)^2 -> expr^2
            => Self::Sq { value: expr },

            // TODO: (expr ^ expr2) ^ expr3 -> expr ^ (expr2 * expr3)

            // TODO: similar to 1^x, x*0 ?

            // TODO:
            //   expr - expr -> 0
            //   expr + expr -> 2 * expr
            //   expr * expr -> expr ^ 2
            //   expr / expr -> 1

            // TODO: ((x + a) + b) -> (x + a), etc

            // TODO: (x - (-b))

            // else do nothing:
            f => f
        };
        if DEBUG { println!("simplify::end: {:?}", new_f) }
        new_f
    }

}



impl ToString for Function {
    fn to_string(&self) -> String {
        match self {
            Self::X => String::from("x"),
            Self::Const { value } => format!("{}", value),
            Self::Zero => String::from("0"),
            Self::One  => String::from("1"),
            Self::Param { name } => format!("{}", name),

            Self::Neg { value } => format!("-{}", value.to_string()),

            Self::Exp { value } => format!("exp({})", value.to_string()),
            Self::Ln  { value } => format!("ln({})", value.to_string()),
            Self::Sqrt{ value } => format!("\\sqrt{{{}}}", value.to_string()),
            Self::Sq  { value } => format!("({})^2", value.to_string()),

            Self::Sin { value } => format!("sin({})", value.to_string()),
            Self::Cos { value } => format!("cos({})", value.to_string()),
            Self::Tan { value } => format!("tan({})", value.to_string()),

            Self::Add { lhs, rhs } => format!("({} + {})", lhs.to_string(), rhs.to_string()),
            Self::Sub { lhs, rhs } => format!("({} - {})", lhs.to_string(), rhs.to_string()),
            Self::Mul { lhs, rhs } => format!("({} * {})", lhs.to_string(), rhs.to_string()),
            Self::Div { lhs, rhs } => format!("({} / {})", lhs.to_string(), rhs.to_string()),
            Self::Pow { lhs, rhs } => format!("({})^{{{}}}", lhs.to_string(), rhs.to_string()),

            Self::Polynomial { degree } => (0..=*degree)
                .map(|i| {
                    format!(
                        "{}{}",
                        PARAMETER_NAMES[i],
                        match i {
                            0 => String::new(),
                            1 => String::from("x"),
                            _ => format!("x^{}", i),
                        }
                    )
                })
                .reduce(|acc, el| format!("{} + {}", acc, el))
                .unwrap(),
            Self::BtrPolynomial { degree } => (0..=*degree)
                .map(|i| {
                    format!(
                        "{}{}/{}!",
                        PARAMETER_NAMES[i],
                        match i {
                            0 => String::new(),
                            1 => String::from("x"),
                            _ => format!("x^{}", i),
                        },
                        i
                    )
                })
                .reduce(|acc, el| format!("{} + {}", acc, el))
                .unwrap(),
        }
    }
}





#[cfg(test)]
mod tests {
    use super::*;

    mod simplify {
        use super::*;
        #[test]
        fn neg_neg_x() {
            assert_eq!(
                Function::X,
                Function::Neg {
                    value: box Function::Neg {
                        value: box Function::X
                    }
                }.simplify()
            );
        }
        #[test]
        fn _4neg_x() {
            assert_eq!(
                Function::X,
                Function::Neg {
                    value: box Function::Neg {
                        value: box Function::Neg {
                            value: box Function::Neg {
                                value: box Function::X
                            }
                        }
                    }
                }.simplify()
            );
        }
        #[test]
        fn _6neg_x() {
            assert_eq!(
                Function::X,
                Function::Neg {
                    value: box Function::Neg {
                        value: box Function::Neg {
                            value: box Function::Neg {
                                value: box Function::Neg {
                                    value: box Function::Neg {
                                        value: box Function::X
                                    }
                                }
                            }
                        }
                    }
                }.simplify()
            );
        }
        #[test]
        fn exp_neg_neg_x() {
            assert_eq!(
                Function::Exp {
                    value: box Function::X
                },
                Function::Exp {
                    value: box Function::Neg {
                        value: box Function::Neg {
                            value: box Function::X
                        }
                    }
                }.simplify()
            );
        }
        #[test]
        fn exp_exp_0() {
            assert_eq!(
                Function::Const {
                    value: 2.718281828459045,
                },
                Function::Exp {
                    value: box Function::Exp {
                        value: box Function::Const { value: 0. }
                    }
                }.simplify()
            );
        }
        #[ignore]
        #[test]
        fn x_x() {
            assert_eq!(
                Function::Sq { value: box Function::X },
                Function::Mul {
                    lhs: box Function::X,
                    rhs: box Function::X
                }.simplify()
            );
        }
    }

    mod from_str {
        use super::*;
        // "literals":
        #[test]
        fn x() {
            assert_eq!(
                Ok(Function::X),
                Function::from_str("x")
            );
        }
        #[test]
        fn const_() {
            assert_eq!(
                Ok(Function::Const { value: 3.14 }),
                Function::from_str("3.14")
            );
        }
        #[test]
        fn zero() {
            assert_eq!(
                Ok(Function::Zero),
                Function::from_str("0")
            );
        }
        #[test]
        fn one() {
            assert_eq!(
                Ok(Function::One),
                Function::from_str("1")
            );
        }
        #[test]
        fn param() {
            for letter in PARAMETER_NAMES {
                assert_eq!(
                    Ok(Function::Param { name: letter }),
                    Function::from_str(&letter.to_string())
                );
            }
        }
        // functions:
        #[test]
        fn neg() {
            assert_eq!(
                Ok(Function::Neg { value: box Function::X }),
                Function::from_str("-x")
            );
        }
        #[test]
        fn exp() {
            assert_eq!(
                Ok(Function::Exp { value: box Function::X }),
                Function::from_str("exp(x)")
            );
        }
        #[test]
        fn ln() {
            assert_eq!(
                Ok(Function::Ln { value: box Function::X }),
                Function::from_str("ln(x)")
            );
        }
        #[test]
        fn sqrt() {
            assert_eq!(
                Ok(Function::Sqrt { value: box Function::X }),
                Function::from_str("sqrt(x)")
            );
        }
        #[test]
        fn sq() {
            assert_eq!(
                Ok(Function::Sq { value: box Function::X }),
                Function::from_str("x^2")
            );
        }
        #[test]
        fn sq_with_brackets() {
            assert_eq!(
                Ok(Function::Sq { value: box Function::X }),
                Function::from_str("(x)^2")
            );
        }
        #[test]
        fn sin() {
            assert_eq!(
                Ok(Function::Sin { value: box Function::X }),
                Function::from_str("sin(x)")
            );
        }
        #[test]
        fn cos() {
            assert_eq!(
                Ok(Function::Cos { value: box Function::X }),
                Function::from_str("cos(x)")
            );
        }
        #[test]
        fn tan() {
            assert_eq!(
                Ok(Function::Tan { value: box Function::X }),
                Function::from_str("tan(x)")
            );
        }
        // operations:
        #[test]
        fn add() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::X,
                    rhs: box Function::Const { value: 3.14 }
                }),
                Function::from_str("x + 3.14")
            );
        }
        #[test]
        fn sub() {
            assert_eq!(
                Ok(Function::Sub {
                    lhs: box Function::X,
                    rhs: box Function::Const { value: 3.14 }
                }),
                Function::from_str("x - 3.14")
            );
        }
        #[test]
        fn mul() {
            assert_eq!(
                Ok(Function::Mul {
                    lhs: box Function::X,
                    rhs: box Function::Const { value: 3.14 }
                }),
                Function::from_str("x * 3.14")
            );
        }
        #[test]
        fn div() {
            assert_eq!(
                Ok(Function::Div {
                    lhs: box Function::X,
                    rhs: box Function::Const { value: 3.14 }
                }),
                Function::from_str("x / 3.14")
            );
        }
        #[test]
        fn pow() {
            assert_eq!(
                Ok(Function::Pow {
                    lhs: box Function::X,
                    rhs: box Function::Const { value: 3.14 }
                }),
                Function::from_str("x ^ 3.14")
            );
        }
        // operations_order
        #[test]
        fn operations_neg_x_sq() {
            assert_eq!(
                Ok(Function::Neg {
                    value: box Function::Sq {
                        value: box Function::X
                    }
                }),
                Function::from_str("-x^2")
            );
        }
        #[test]
        fn operations_neg_lb_x_rb_sq() {
            assert_eq!(
                Ok(Function::Neg {
                    value: box Function::Sq {
                        value: box Function::X
                    }
                }),
                Function::from_str("-(x)^2")
            );
        }
        #[test]
        fn operations_lb_neg_x_rb_sq() {
            assert_eq!(
                Ok(Function::Sq {
                    value: box Function::Neg {
                        value: box Function::X
                    }
                }),
                Function::from_str("(-x)^2")
            );
        }
        #[test]
        fn operations_neg_a_add_b() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::Neg {
                        value: box Function::Param { name: 'a' }
                    },
                    rhs: box Function::Param { name: 'b' }
                }),
                Function::from_str("-a+b")
            );
        }
        #[test]
        fn operations_neg_sin_a_add_cos_b() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::Neg {
                        value: box Function::Sin {
                            value: box Function::Param { name: 'a' }
                        }
                    },
                    rhs: box Function::Cos {
                        value: box Function::Param { name: 'b' }
                    }
                }),
                Function::from_str("-sin(a)+cos(b)")
            );
        }
        #[test]
        fn operations_neg_a_mul_b() {
            assert_eq!(
                Ok(Function::Neg {
                    value: box Function::Mul {
                        lhs: box Function::Param { name: 'a' },
                        rhs: box Function::Param { name: 'b' }
                    }
                }),
                Function::from_str("-a*b")
            );
        }
        #[test]
        fn operations_order_add_mul() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::Const { value: 145. },
                    rhs: box Function::Mul {
                        lhs: box Function::Const { value: 42. },
                        rhs: box Function::Const { value: 3.14 }
                    }
                }),
                Function::from_str("145 + 42 * 3.14")
            );
        }
        #[test]
        fn operations_order_mul_add() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::Mul {
                        lhs: box Function::Const { value: 145. },
                        rhs: box Function::Const { value: 42. }
                    },
                    rhs: box Function::Const { value: 3.14 }
                }),
                Function::from_str("145 * 42 + 3.14")
            );
        }
        #[test]
        fn operations_order_lb_add_rb_mul() {
            assert_eq!(
                Ok(Function::Mul {
                    lhs: box Function::Add {
                        lhs: box Function::Const { value: 145. },
                        rhs: box Function::Const { value: 42. }
                    },
                    rhs: box Function::Const { value: 3.14 }
                }),
                Function::from_str("(145 + 42) * 3.14")
            );
        }
        #[test]
        fn operations_order_add_lb_mul_rb() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::Const { value: 145. },
                    rhs: box Function::Mul {
                        lhs: box Function::Const { value: 42. },
                        rhs: box Function::Const { value: 3.14 }
                    }
                }),
                Function::from_str("145 + (42 * 3.14)")
            );
        }
        #[test]
        fn operations_order_lb_mul_rb_add() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::Mul {
                        lhs: box Function::Const { value: 145. },
                        rhs: box Function::Const { value: 42. }
                    },
                    rhs: box Function::Const { value: 3.14 }
                }),
                Function::from_str("(145 * 42) + 3.14")
            );
        }
        #[test]
        fn operations_order_mul_lb_add_rb() {
            assert_eq!(
                Ok(Function::Mul {
                    lhs: box Function::Const { value: 145. },
                    rhs: box Function::Add {
                        lhs: box Function::Const { value: 42. },
                        rhs: box Function::Const { value: 3.14 }
                    }
                }),
                Function::from_str("145 * (42 + 3.14)")
            );
        }
        #[test]
        fn operations_order_sub_sub() {
            assert_eq!(
                Ok(Function::Sub {
                    lhs: box Function::Sub {
                        lhs: box Function::Const { value: 145. },
                        rhs: box Function::Const { value: 42. }
                    },
                    rhs: box Function::Const { value: 3.14 }
                }),
                Function::from_str("145 - 42 - 3.14")
            );
        }
        #[test]
        fn operations_order_add_mul_pow() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::Const { value: 145. },
                    rhs: box Function::Mul {
                        lhs: box Function::Const { value: 42. },
                        rhs: box Function::Pow {
                            lhs: box Function::Const { value: 3.14 },
                            rhs: box Function::Const { value: 2.71 }
                        }
                    },
                }),
                Function::from_str("145 + 42 * 3.14 ^ 2.71")
            );
        }
        #[test]
        fn operations_order_add_pow_mul() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::Const { value: 145. },
                    rhs: box Function::Mul {
                        lhs: box Function::Pow {
                            lhs: box Function::Const { value: 42. },
                            rhs: box Function::Const { value: 3.14 }
                        },
                        rhs: box Function::Const { value: 2.71 }
                    },
                }),
                Function::from_str("145 + 42 ^ 3.14 * 2.71")
            );
        }
        #[test]
        fn operations_order_mul_add_pow() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::Mul {
                        lhs: box Function::Const { value: 145. },
                        rhs: box Function::Const { value: 42. }
                    },
                    rhs: box Function::Pow {
                        lhs: box Function::Const { value: 3.14 },
                        rhs: box Function::Const { value: 2.71 }
                    }
                }),
                Function::from_str("145 * 42 + 3.14 ^ 2.71")
            );
        }
        #[test]
        fn operations_order_mul_pow_add() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::Mul {
                        lhs: box Function::Const { value: 145. },
                        rhs: box Function::Pow {
                            lhs: box Function::Const { value: 42. },
                            rhs: box Function::Const { value: 3.14 }
                        }
                    },
                    rhs: box Function::Const { value: 2.71 },
                }),
                Function::from_str("145 * 42 ^ 3.14 + 2.71")
            );
        }
        #[test]
        fn operations_order_pow_add_mul() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::Pow {
                        lhs: box Function::Const { value: 145. },
                        rhs: box Function::Const { value: 42. }
                    },
                    rhs: box Function::Mul {
                        lhs: box Function::Const { value: 3.14 },
                        rhs: box Function::Const { value: 2.71 }
                    }
                }),
                Function::from_str("145 ^ 42 + 3.14 * 2.71")
            );
        }
        #[test]
        fn operations_order_pow_mul_add() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::Mul {
                        lhs: box Function::Pow {
                            lhs: box Function::Const { value: 145. },
                            rhs: box Function::Const { value: 42. }
                        },
                        rhs: box Function::Const { value: 3.14 }
                    },
                    rhs: box Function::Const { value: 2.71 },
                }),
                Function::from_str("145 ^ 42 * 3.14 + 2.71")
            );
        }
        // complex:
        #[test]
        fn from_str() {
            assert_eq!(
                Ok(Function::Neg {
                    value: box Function::Sin {
                        value: box Function::X
                    }
                }),
                Function::from_str("-sin(x)")
            );
        }
        #[test]
        fn complex_two_gausses() {
            assert_eq!(
                Ok(Function::Add {
                    lhs: box Function::Param { name: 'h' },
                    rhs: box Function::Add {
                        lhs: box Function::Mul {
                            lhs: box Function::Param { name: 'a' },
                            rhs: box Function::Exp {
                                value: box Function::Neg {
                                    value: box Function::Sq {
                                        value: box Function::Div {
                                            lhs: box Function::Sub {
                                                lhs: box Function::X,
                                                rhs: box Function::Param { name: 'm' }
                                            },
                                            rhs: box Function::Param { name: 's' }
                                        }
                                    }
                                }
                            }
                        },
                        rhs: box Function::Mul {
                            lhs: box Function::Param { name: 'b' },
                            rhs: box Function::Exp {
                                value: box Function::Neg {
                                    value: box Function::Sq {
                                        value: box Function::Div {
                                            lhs: box Function::Sub {
                                                lhs: box Function::X,
                                                rhs: box Function::Param { name: 'n' }
                                            },
                                            rhs: box Function::Param { name: 't' }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }),
                Function::from_str("h + a*exp(-((x-m)/s)^2) + b*exp(-((x-n)/t)^2)")
            );
        }
        #[test]
        fn complex() {
            assert_eq!(
                Ok(Function::Mul {
                    lhs: box Function::Mul {
                        lhs: box Function::One,
                        rhs: box Function::X
                    },
                    rhs: box Function::Add {
                        lhs: box Function::Div {
                            lhs: box Function::Param { name: 'c' },
                            rhs: box Function::Add {
                                lhs: box Function::Neg {
                                    value: box Function::Mul {
                                        lhs: box Function::X,
                                        rhs: box Function::Param { name: 't' }
                                    }
                                },
                                rhs: box Function::Add {
                                    lhs: box Function::Param { name: 'u' },
                                    rhs: box Function::Sin {
                                        value: box Function::Mul {
                                            lhs: box Function::Param { name: 'v' },
                                            rhs: box Function::X
                                        }
                                    }
                                }
                            }
                        },
                        rhs: box Function::Param { name: 'm' }
                    }
                }),
                Function::from_str("(1*x) * ((c/(-(x*t)+(u+sin(v*x))))+m)")
            );
        }
        #[test]
        fn complex_2() {
            assert_eq!(
                Ok(Function::Mul {
                    lhs: box Function::Pow {
                        lhs: box Function::Pow {
                            lhs: box Function::Div {
                                lhs: box Function::Exp {
                                    value: box Function::X
                                },
                                rhs: box Function::X,
                            },
                            rhs: box Function::Param { name: 'w' }
                        },
                        rhs: box Function::Param { name: 'q' }
                    },
                    rhs: box Function::Mul {
                        lhs: box Function::X,
                        rhs: box Function::Param { name: 'v' }
                    }
                }),
                Function::from_str("((exp(x) / x)^(w))^(q) * (x * v)")
            );
        }
    }

    mod eval {
        use super::*;
        #[test]
        fn x() {
            assert_eq!(
                2.5,
                Function::from_str("x").unwrap()
                    .eval(2.5, &Params::new()),
            );
        }
        #[test]
        fn complex() {
            assert_eq!(
                11.479758672104968,
                Function::from_str("x + 2*a*(x+1)^2 - sin(x+1) + exp(3*x)/exp(x)").unwrap()
                    .eval(1., &Params::from_array([('a', 0.5)])),
            );
        }
    }
}

