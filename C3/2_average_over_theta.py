#average_over_theta.py: 計算したstopping or inteference_functionをある区間thetasで区切り、その平均を取る(row: charge, col: theta)
#n=3
import os
import sys

#平均を取る区間
thetas = [0, 30, 60, 90]

working_dir =  r'results'

target = 'Gly'
input_dirname = working_dir
input_filename = r'E=900keV_atom_C3_linear_{}.txt'.format(target)
output_dirname = "averaged"

def main():
    #二次元配列として読み込み
    with open(os.path.join(input_dirname, input_filename)) as f:
        dat = [v.split() for v in f.read().split('\n') if len(v) != 0]

    #output_path
    output_dirpath = os.path.join(input_dirname, output_dirname)
    os.makedirs(output_dirpath, exist_ok=True)
    output_path  = os.path.join(output_dirpath, input_filename[:-4] + '_ave_para={}deg-.txt'.format(thetas[0]))
    output_ratio_path  = os.path.join(output_dirpath, input_filename[:-4] + '_ave_ratio.txt')

    output_data = []
    header = []
    for i in range(len(thetas)-1):
        header.append(str((thetas[i] + thetas[i+1])/2))
    output_data.append(header)
    
    output_ratio_data = []

    for row in range(1, len(dat)):
        #0列目:charge pair
        output_line = []
        output_ratio_line = []

        #thetaの区間数についてループ
        for i in range(len(thetas)-1):
            #区間の定義
            t_start = thetas[i]
            t_end = thetas[i+1]
            #出力の一行
            tmp = []
            #thetaについてループして、区間に含まれるデータを列挙する
            for col in range(1, len(dat[0])):
                #0行目のthetaを参照
                if t_start <= int(dat[0][col]) < t_end:
                    tmp.append(float(dat[row][col]))
            #平均を取り、append
            ave = sum(tmp)/len(tmp)
            output_line.append(str(ave))

        #出力行を追加
        output_data.append(output_line)
        output_ratio_line.append(str(float(output_line[2])/float(output_line[0])))
        output_ratio_data.append(output_ratio_line)
    
    #保存
    with open(output_path, 'w') as f:
        #write
        for i in range(len(output_data)):
            f.write('\t'.join(output_data[i])+'\n')
    with open(output_ratio_path, 'w') as f:
        #write
        for i in range(len(output_ratio_data)):
            f.write('\t'.join(output_ratio_data[i])+'\n')
            
    print(input_filename + '---done')
    print('finished!')

if __name__ == '__main__':
    main()