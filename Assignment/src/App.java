import java.math.BigInteger;
import java.util.Random;

public class App {

    /* Ham GCD su dung phuong phap Euclid de tim uoc so chung lon nhat cua a va b */
    static BigInteger GCD(BigInteger a, BigInteger b) {
        while (true) {
            if (a.compareTo(BigInteger.ZERO) == 0) return b;
            if (b.compareTo(BigInteger.ZERO) == 0) return a;
            if (a.compareTo(b) == 1) {
                a = a.mod(b);
            } else {
                b = b.mod(a);
            }
        }
    }

    /* Ham Inverse su dung de quy de tim gia tri x va y phu hop voi dieu kien da cho:
        - Neu b bang 1, ta tra ve gia tri x bang 1 va y bang a - 1. Dieu nay dua tren tinh chat: a1 - b(a-1) = a - b*a + b = 1.
        - Neu a bang 1, ta tra ve gia tri x bang 1 va y bang 0. Dieu nay dua tren tinh chat: a1 - b0 = a.
        - Neu a lon hon b, ta thuc hien phep toan % de lay phan du cua a chia cho b, va tiep tuc goi de quy ham Inverse voi cac tham so a la phan du tren va b la b. 
        Sau đo, ta tinh toan gia tri x va y nhu sau:
            x bang gia tri x duoc tra ve boi ham de quy.
            y bang a/b*x + y, trong do a/b la phep chia nguyen cua a cho b, va y la gia tri y duoc tra ve boi ham de quy.
        Neu a nho hon b, ta thuc hien phep toan % de lay phan du cua b chia cho a, va tiep tuc goi de quy ham Inverse voi cac tham so a la a va b la phan du tren. Sau đo, ta tinh toan gia tri x va y nhu sau:
            x bang b/a*y + x, trong do b/a la phep chia nguyen cua b cho a, va x la gia tri x duoc tra ve boi ham de quy.
            y bang gia tri y duoc tra ve boi ham de quy.
        Cuoi cung, ham Inverse tra ve mot mang BigInteger chua gia tri x va y tuong ung.
    */
    static BigInteger[] Inverse(BigInteger a, BigInteger b) { 
        if (b.compareTo(BigInteger.ONE) == 0) {
            return new BigInteger[] { 
                BigInteger.ONE, 
                a.subtract(BigInteger.ONE) 
            };
        }
        if (a.compareTo(BigInteger.ONE) == 0) {
            return new BigInteger[] {
                BigInteger.ONE, 
                BigInteger.ZERO
            };
        }
        if (a.compareTo(b) != -1) { 
            BigInteger[] vals = Inverse(a.mod(b), b);
            return new BigInteger[] {
                vals[0], 
                a.divide(b).multiply(vals[0]).add(vals[1])
            };
        }  
        BigInteger[] vals = Inverse(a, b.mod(a));
        return new BigInteger[] {
            b.divide(a).multiply(vals[1]).add(vals[0]),
            vals[1]
        };
    }

    /* Thuc hien tinh a^b mod n, voi a, b, n la cac so nguyen duong duoc bieu dien bang doi tuong BigInteger.
        Cu the, ham duoc thuc hien nhu sau:
            - Neu b = 0, tra ve gia tri 1.
            - Goi ham de quy de tinh a^(b/2) mod n, luu ket qua vao bien res.
            - Neu b chan, tra ve gia tri res^2 mod n.
            - Neu b le, tra ve gia tri res^2 * a mod n.
    */
    static BigInteger Pow(BigInteger a, BigInteger b, BigInteger n) {
        BigInteger two = new BigInteger("2");
        if (b.compareTo(BigInteger.ZERO) == 0) {
            return BigInteger.ONE;
        }
        BigInteger res = Pow(a, b.divide(two), n);
        if (b.mod(two).compareTo(BigInteger.ZERO) == 0) {
            return res.multiply(res).mod(n);
        }
        return res.multiply(res).mod(n).multiply(a).mod(n);
    }

    /*  Tim gia tri d sao cho ed = 1 (mod phi(n)), trong do e la mot so nguyen to cung nhau voi phi(n), 
        p va q la hai so nguyen to khac nhau, n = pq, va phi(n) = (p-1)*(q-1).
            - Dau vao: e, p, q la ba tham so kieu BigInteger.
            - Dau ra: mot gia tri d kieu BigInteger thoa man ed = 1 (mod phi(n)).
        De tim gia tri d, ta su dung ham Inverse() de tim gia tri nguoc mod cua e theo phi(n), nghia la tim cap so nguyen x, y sao cho ex - y*phi(n) = 1. 
        Ham Inverse() duoc dinh nghia truoc do trong doan code. Sau do, ta chi can lay gia tri x cua cap so tra ve tu ham Inverse() la gia tri cua d.
    */
    static BigInteger PhiInv(BigInteger e, BigInteger p, BigInteger q) {
        BigInteger pminus1 = p.subtract(BigInteger.ONE);
        BigInteger qminus1 = q.subtract(BigInteger.ONE);
        BigInteger phi = pminus1.multiply(qminus1);
        BigInteger[] arr = Inverse(e, phi);
        BigInteger d = arr[0];
        return d;
    }

    /*  Sinh mot so nguyen ngau nhien nho hon mot so nguyen duong n cho truoc. 
        No su dung doi tuong cua lop java.util.Random de sinh so ngau nhien va tra ve mot doi tuong BigInteger voi gia tri nho hon n. 
    */
    static BigInteger RandomBigInteger(BigInteger m, BigInteger n) { //Tao bo khoa ngau nhien 2 so nguyen to lon
        Random randNum = new Random();
        int len = n.bitLength();
        BigInteger res = new BigInteger(len, randNum);
        while ((res.compareTo(n) != -1) || (res.compareTo(m) == -1)) {
            randNum = new Random();
            res = new BigInteger(len, randNum);
        }
        return res;
    }

    /* Thuat toan kiem tra tinh nguyen to cua so nguyen duong n bang phuong phap Miller Rabin. 
        Dau vao cua ham la so nguyen duong n can kiem tra tinh nguyen to va mot so nguyen duong a duoc goi la chung nhan. 
        Neu ham tra ve gia tri true thi n co the la so nguyen to, con neu tra ve gia tri false thi n chac chan khong phai la so nguyen to.

        Thuat toan bat dau bang viec tinh gia tri cua s, la mot so nguyen duong sao cho (n-1) co the viet duoi dang 2^s * d. 
        Sau do, ta tinh gia tri cua d. Tiep theo, ta tinh gia tri cua x bang cach tinh a^d mod n. Ta lap lai buoc nay cho den khi i bang s-1. 
        Sau moi lan lap, ta tinh gia tri cua y bang x^2 mod n. 
        Neu tai bat ky lan lap nao, y = 1 va x khac 1 va n-1, thi ta ket luan rang n khong phai la so nguyen to va tra ve gia tri false. 
        Neu y khac 1 o cuoi cung, ta cung ket luan rang n khong phai la so nguyen to va tra ve gia tri false. 
        Neu khong co truong hop tren xay ra, ta ket luan rang n co the la so nguyen to va tra ve gia tri true.
    */
    static boolean Miller_Rabin (BigInteger n, BigInteger a) {
        BigInteger nminus1 = n.subtract(BigInteger.ONE);
        BigInteger two = new BigInteger("2");
        BigInteger copy = n.subtract(BigInteger.ONE);
        int s = 0;
        while (copy.mod(two).compareTo(BigInteger.ZERO) == 0) {
            s++;
            copy = copy.divide(two);
        }
        BigInteger y = BigInteger.ONE;
        BigInteger d = nminus1.divide(two.pow(s));
        BigInteger x = Pow(a, d, n);
        for (int i = 0; i < s; i++) {
            y = Pow(x, two, n);
            if ((y.compareTo(BigInteger.ONE) == 0) && (x.compareTo(BigInteger.ONE) != 0) && (x.compareTo(nminus1) != 0)){
                return false;
            }
            x = y;
        }
        if (y.compareTo(BigInteger.ONE) != 0) {
            return false;
        }
        return true;
    }

    /*  Kiem tra so nguyen to Miller-Rabin de kiem tra tinh nguyen to cua mot so nguyen lon n. 
        Ham nay kiem tra n bang cach chon ngau nhien 15 so nguyen to dau tien va kiem tra n co chia het cho bat ky so nguyen to do khong. 
        Sau do, ham chon ngau nhien 15 so a va kiem tra n bang cach su dung phuong phap kiem tra so nguyen to Miller-Rabin voi so chung nhan a da chon. 
        Neu n khong vuot qua bat ky buoc kiem tra nao, ham se tra ve false, cho biet rang n la mot so hop so. 
        Neu n vuot qua tat ca cac buoc kiem tra, ham se tra ve true, cho biet rang n co xac suat cao la mot so nguyen to.

        Ham cung su dung ham RandomBigInteger de tao so nguyen to ngau nhien, va ham GCD de kiem tra xem a co la uoc chung lon nhat voi n hay khong. 
        Ngoai ra, ham con su dung danh sach cac so nguyen to dau tien de kiem tra xem n co chia het cho bat ky so nguyen to nao trong danh sach do khong, 
        va tra ve false neu n chia het cho mot trong so do.
    */
    static boolean isPrime(BigInteger n) {
        boolean bool = true;
        int t = 15;
        BigInteger a = new BigInteger("1");
        BigInteger one = a;

        for (int i = 0; i < t; i++) {
            a = RandomBigInteger(one, n);
            if (GCD(a, n).compareTo(BigInteger.ONE) != 0) {
                return false;
            }
            bool = Miller_Rabin(n, a);
            if (bool == false) {
                return bool;
            }
        }

        String[] prime = { "2", "3", "5", "7", "11", "13", "17", "19", "23", "29", "31", "37", "41", "43", "47", "53", 
            "59", "61", "67", "71", "73", "79", "83", "91", "101", "103", "107", "109", "113", "127", "131", "137", "139", "149", "151", "157", "163" };
        for (int i = 0; i < 37; i++) {
            a = new BigInteger(prime[i]);
            if (n.mod(a).compareTo(BigInteger.ZERO) == 0) {
                return false;
            }
        }
        return bool;
    }

    /*  Tao ra mot so nguyen to ngau nhien voi do dai bit duoc chi dinh:
        - Dau tien, phuong thuc su dung lop Random de tao mot so nguyen ngau nhien voi do dai bit chi dinh va luu ket qua vao bien res.
        - Sau do, phuong thuc kiem tra xem so nguyen res co chia het cho 2 hay khong. Neu co, no tang res len 1 de dam bao rang res la so le.
        - Tiep theo, phuong thuc kiem tra xem res co phai la so nguyen to hay khong bang cach goi phuong thuc isPrime. 
        Neu res khong phai la so nguyen to, phuong thuc se thuc hien buoc tiep theo.
        - Neu res khong phai la so nguyen to, phuong thuc se kiem tra xem res co chia het cho 6 hay khong. 
        Neu res = 1 (mod 6), nghia la res co the co dang 6k + 1, vi vay phuong thuc tang res len 4 de no co dang 6k + 5. 
        Neu res khong chia het cho 6, phuong thuc tang res len 2 de dam bao rang res la so le.
        - Sau khi tang gia tri cua res, phuong thuc kiem tra lai xem res co phai la so nguyen to hay khong bang cach goi phuong thuc isPrime. 
        Phuong thuc lap lai cac buoc 4 va 5 cho den khi tim duoc mot so nguyen to. Khi do, phuong thuc tra ve so nguyen to do la ket qua.
        - Phuong thuc su dung kiem tra Miller-Rabin voi so luong lap lai la 15 lan de kiem tra tinh nguyen to cua mot so n. 
        Neu n khong phai la so nguyen to, kha nang de phuong thuc tra ve ket qua sai la rat nho, nho hon 1/10^6.
    */
    static BigInteger getPrime(int bit) {
        Random randNum = new Random();
        BigInteger res = new BigInteger(bit, randNum);
        BigInteger two = new BigInteger("2");
        BigInteger four = new BigInteger("4");
        BigInteger six = new BigInteger("6");

        boolean t = isPrime(res);
        if (res.mod(two).compareTo(BigInteger.ZERO) == 0) {
            res = res.add(BigInteger.ONE);
        }
        while (t == false) {
            if (res.mod(six).compareTo(BigInteger.ONE) == 0) {
                res = res.add(four);
            } else {
                res = res.add(two);
            }
            t = isPrime(res);
        }
        return res;
    }

    /*  Tao ra ba so nguyen to manh (strong prime) p, q, va a, duoc su dung trong he mat RSA. 
        Cac so nguyen to manh duoc dinh nghia la cac so nguyen to p va q sao cho (p-1)/2 va (q-1)/2 deu la so nguyen to. So a duoc su dung de tao khoa bi mat.

        Dau tien, chung ta khoi tao cac bien can thiet, bao gom cac doi tuong BigInteger va mot so hang so. 
        Sau do, chung ta tao ra mot so ngau nhien res trong khoang giua lo va hi.

        Sau do, chung ta tim p0 va p1 bang cach su dung ham getPrime de tao ra hai so nguyen to co do dai 128 bit. 
        Chung ta xac dinh rang e khong chia het cho p0 va p1 phai khac nhau.

        Tiep theo, chung ta tinh toan gia tri crt, la gia tri thoa man crt = 1 (mod p0) va crt = -1 (mod p1). 
        Chung ta tinh toan hai gia tri crt1 va crt2 tuong ung voi p0 va p1, sau do cong chung de co duoc gia tri crt. 
        Chung ta sau do tinh toan gia tri res sao cho res = crt (mod p0*p1), res-1 chia het cho p0 va res+1 chia het cho p1.

        Cuoi cung, chung ta kiem tra tung so nguyen to co the co bang cach thu cong them vao gia tri res mot luong tang dan cua increment, voi increment bang p0p12. 
        Neu so nguyen to duoc tim thay, chung ta tra ve cac so nguyen to manh p, q va a duoc tinh toan tu cac gia tri res, p0 va p1.
    */
    static BigInteger[] getstrongPrime(int bit, BigInteger e) {
        BigInteger two = new BigInteger("2");
        BigInteger lo = new BigInteger( "94807519081091767274915958072241511435209662959" +
            "194056066952784574407133867106560704159305496966045119348817003874896742069226" +
            "49072631363246958077487951940"
        );
        BigInteger hi = two.pow(512).subtract(BigInteger.ONE);
        BigInteger res = RandomBigInteger(lo, hi);
        boolean result = false;
        BigInteger p0 = getPrime(128);
        while (e.mod(p0).compareTo(BigInteger.ZERO) == 0) {
            p0 = getPrime(128);
        }
        BigInteger p1 = getPrime(128);
        while ((p1.compareTo(p0) == 0) || (e.mod(p1).compareTo(BigInteger.ZERO) == 0)) {
            p1 = getPrime(128);
        }
        p1 = p1.multiply(two);

        BigInteger increment = p0.multiply(p1);

        BigInteger[] inv1 = Inverse(p1, p0);
        BigInteger crt1 = inv1[0]; 
        crt1 = crt1.multiply(p1);

        BigInteger[] inv2 = Inverse(p0, p1);
        BigInteger crt2 = p1.subtract(inv2[0]); 
        crt2 = crt2.multiply(p0);

        BigInteger crt = crt1.add(crt2).add(increment);
        BigInteger resmod = res.mod(increment);
        res = res.add(crt.subtract(resmod).mod(increment));
        increment = increment.multiply(two);

        while (true) {
            boolean possible_prime = true;
            if (possible_prime) {
                if (GCD(e, res.subtract(BigInteger.ONE)).compareTo(BigInteger.ONE) != 0) { 
                    possible_prime=false;
                }
            }
            if (possible_prime) {
                result = isPrime(res);
                if (result == true) {
                    break;
                }
            }
            res = res.add(increment);
        }
        return new BigInteger[] {res, p0, p1.divide(two)};
    }

    static BigInteger Encrypt(BigInteger m, BigInteger e, BigInteger p, BigInteger q) {
        BigInteger N = p.multiply(q);
        return Pow(m, e, N);
    }

    static BigInteger Decrypt(BigInteger c, BigInteger d, BigInteger p, BigInteger q) {
        BigInteger N = p.multiply(q);
        return Pow(c, d, N);
    }

    public static void main(String args[]) throws Exception {
        BigInteger two = new BigInteger("2");
        
        /*  Buoc dau tien la sinh ra khoa cong khai va khoa bi mat. Trong do, khoa cong khai duoc tao bang cach sinh ra mot so nguyen to ngau nhien co do dai la 1024 bit, goi la e. 
            Sau do, kiem tra xem e co chia het cho 2 hay khong, neu co thi tang gia tri cua e len mot don vi de dam bao e la so le.
        */ 
        System.out.println("Sinh ra e");
        Random rand = new Random();
        BigInteger e = new BigInteger(1024, rand);
        if (e.mod(two).compareTo(BigInteger.ZERO) == 0) {
            e = e.add(BigInteger.ONE);
        }
        System.out.println("\t" + e); 
        System.out.println("\n");

        System.out.println("Sinh ra Strong-prime p:");
        BigInteger p = new BigInteger("512");
        BigInteger[] list = getstrongPrime(512, e);
        p = list[0];
        BigInteger p0 = list[1];
        BigInteger p1 = list[2];
        System.out.println("\t" + p);
        System.out.println("p - 1 chia het cho mot Large-prime: " + p0);
        System.out.println("p + 1 chia het cho mot Large-prime: " + p1);
        System.out.println("\n");

        /*  Tao ra hai so nguyen to manh la p va q. Hai so nguyen to nay co do dai la 512 bit va duoc tao ra bang phuong phap chon ngau nhien mot so nguyen to lon roi kiem tra tinh manh cua no 
            bang cach kiem tra tinh chia het cua p-1 va p+1 cho mot so nguyen to lon, tuong tu voi q-1 va q+1.
        */
        System.out.println("Sinh ra Strong-prime q");
        BigInteger q = new BigInteger("512");
        BigInteger[] list2 = getstrongPrime(512, e);
        q = list2[0];
        BigInteger q0 = list2[1];
        BigInteger q1 = list2[2];
        System.out.println("\t" + q);
        System.out.println("q - 1 chia het cho mot Large-prime: " + q0);
        System.out.println("q + 1 chia het cho mot Large-prime: " + q1);
        System.out.println("\n");

        System.out.println("Message:"); //Sinh ra noi dung nguoi gui muon gui di voi dieu kien la do dai cua noi dung nho hon do dai cua p*q
        Random randNum = new Random();
        BigInteger m = new BigInteger(1024, randNum);
        //Ensure that m < p*q
        m = m.mod(p.multiply(q));
        System.out.println(m); 
        System.out.println("\n");

        /*  Ma hoa Encrypt tin nhan bang cach su dung khoa cong khai e, va hai so nguyen to p va q da duoc tao ra truoc do.*/
        System.out.println("Encryption:");
        BigInteger c = Encrypt(m, e, p, q);
        System.out.println(c);  
        System.out.println("\n");

        /*  Doan code tinh toan khoa bi mat d bang cach tinh nghich dao modulo cua e theo cong thuc d = e^(-1) mod phi(n), trong do phi(n) duoc tinh bang (p-1)*(q-1).
            Su dung ham PhiInv ben tren.
         */
        System.out.println("Inverse modulo phi:");
        BigInteger d = PhiInv(e, p, q);
        System.out.println(d);  
        System.out.println("\n");

        /*  Doan code giai ma tin nhan bang cach su dung khoa bi mat d va hai so nguyen to p va q da duoc tao ra truoc do. Ket qua giai ma se giong voi tin nhan ban dau.
            Giai ma Decrypt.
        */
        System.out.println("Decryption:");
        m = Decrypt(c, d, p, q);
        assert(m.compareTo(p.multiply(q)) == -1);
        System.out.println(m);
        System.out.println("\n");
    }
}