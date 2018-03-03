class Runge
  require 'narray'
  require 'pp'
  H, X0, Xn, S, Xin, Yin, DXin, DYin = 0.0001, 0.0, 10.0, 5, 0.0, 1.0, 0.0, -1.0
  C = NArray[1.0/4.0, 3.0/4.0, 11.0/20.0, 1.0/2.0, 1.0]
  #B = NArray[25.0/24.0, -49.0/48.0, 125.0/16.0, -85.0/12.0, 1.0/4.0]
  B = NArray[59.0/48.0, -17.0/96.0, 225.0/32.0, -85.0/12.0, 0.0]
  A = NArray[[1.0/4.0, 0.0, 0.0, 0.0, 0.0],
             [1.0/2.0, 1.0/4.0, 0.0, 0.0, 0.0],
             [17.0/50.0, -1.0/25.0, 1.0/4.0, 0.0, 0.0],
             [371.0/1360.0, -137.0/2720.0, 15.0/544.0, 1.0/4.0, 0.0],
             [25.0/24.0, -49.0/48.0, 125.0/16.0, -85.0/12.0, 1.0/4.0]]
  def self.f1(u, v)
    v
  end
  def self.summ(u, v, p)
    su = 0.0
    S.times{|i| su += u[p, i] * v[i]}
    su
  end
  def self.f2(u, v)
    10.0 * v + 11.0 * u
  end
  y = NArray.float(((Xn - X0) / H).to_int + 1)
  y[0] = Yin
  z = NArray.float(((Xn - X0) / H).to_int + 1)
  z[0] = DYin
  #  i from 1 to (Xn - X0) / H
  f = File.new('lib.txt', 'w')
  f.puts y[0]
  (((Xn - X0) / H).to_int).times do |i|
    k1 = NArray.float(S)
    k2 = NArray.float(S)
    S.times do |j|
      su = 0.0
      S.times{|k| su += A[k, j] * k1[k]}
      k1[j] = f1(y[i] + H * su, z[i] + H * su) / (1.0 - H * f1(A[j, j], A[j, j]))
      su = 0.0
      S.times{|k| su += A[k, j] * k2[k]}
      k2[j] = f2(y[i] + H * su, z[i] + H * su) / (1.0 - H * f2(A[j, j], A[j, j]))
    end
    y[i + 1] = y[i]
    z[i + 1] = z[i]
    S.times do |j|
      y[i + 1] += H * B[j] * k1[j]
      z[i + 1] += H * B[j] * k2[j]
    end
    f.puts y[i + 1]
  end
  f.close
end
