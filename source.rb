class Runge
  require 'narray'
  require 'pp'
  H, X0, Xn, S, Xin, Yin, DXin, DYin, C1, C2, C3, C4 = 0.0001, 0.0, 10.0, 5, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 11.0, 10.0, 0.0, 1.0, 11.0, 10.0
  C = NArray[1.0/4.0, 3.0/4.0, 11.0/20.0, 1.0/2.0, 1.0]
  #B = NArray[25.0/24.0, -49.0/48.0, 125.0/16.0, -85.0/12.0, 1.0/4.0]
  B = NArray[59.0/48.0, -17.0/96.0, 225.0/32.0, -85.0/12.0, 0.0]
  A = NArray[[1.0/4.0, 0.0, 0.0, 0.0, 0.0],
             [1.0/2.0, 1.0/4.0, 0.0, 0.0, 0.0],
             [17.0/50.0, -1.0/25.0, 1.0/4.0, 0.0, 0.0],
             [371.0/1360.0, -137.0/2720.0, 15.0/544.0, 1.0/4.0, 0.0],
             [25.0/24.0, -49.0/48.0, 125.0/16.0, -85.0/12.0, 1.0/4.0]]
  def self.f1(u, v)
    C1 * u + C2 * v
  end
  def self.f2(u, v)
    C3 * u + C4 * v
  end
  y, z = NArray.float(((Xn - X0) / H).to_int + 1), NArray.float(((Xn - X0) / H).to_int + 1)
  y[0], z[0] = Yin, DYin
  #  i from 1 to (Xn - X0) / H
  p = File.new('lib.txt', 'w')
  p.puts Xin.to_s + ";" + y[0].to_s
  (((Xn - X0) / H).to_int).times do |i|
    k1, k2 = NArray.float(S), NArray.float(S)
    S.times do |j|
      su = 0.0
      su1 = 0.0
      S.times do |k|
        su += A[k, j] * k1[k]
        su1 += A[k, j] * k2[k]
      end
      f = f1(y[i] + H * su, z[i] + H * su1)
      g = f2(y[i] + H * su, z[i] + H * su1)
      k1[j] = (f * (1 - C4 * H * A[j, j]) + C2 * H * A[j, j] * g) /
              ((1.0 - C4 * H * A[j, j]) * (1.0 - C1 * H * A[j, j]) -
              C2 * C3 * H * H * A[j, j] * A[j, j])
      k2[j] = ((1 - C1 * H * A[j, j]) * k1[j] - f) / (C2 * H * A[j, j])
    end
    y[i + 1], z[i + 1] = y[i], z[i]
    S.times do |j|
      y[i + 1] += H * B[j] * k1[j]
      z[i + 1] += H * B[j] * k2[j]
    end
    p.puts ((i + 1) * H).to_s + ";" + y[i + 1].to_s
  end
  p.close
end
