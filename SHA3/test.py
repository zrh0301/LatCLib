import hashlib
import sys

def run_tests():
    print(f"{'Algorithm':<20} | {'Input':<10} | {'Output (First 16 Bytes/32 Hex chars)':<40}")
    print("-" * 80)

    # 1. 测试 SHAKE128 (空输入)
    shake = hashlib.shake_128()
    # update nothing
    output_hex = shake.hexdigest(16)
    print(f"{'SHAKE128':<20} | {'(empty)':<10} | {output_hex}")

    # 2. 测试 SHAKE256 (输入 "abc")
    shake = hashlib.shake_256()
    shake.update(b"abc")
    output_hex = shake.hexdigest(16)
    print(f"{'SHAKE256':<20} | {'abc':<10} | {output_hex}")

    # 3. 测试 SHAKE128 (输入 "abc")
    # 这对应你在 C main.c 中的第三个测试 (Stream: "a" + "bc")
    shake = hashlib.shake_128()
    shake.update(b"abc")
    output_hex = shake.hexdigest(16)
    print(f"{'SHAKE128':<20} | {'abc':<10} | {output_hex}")

if __name__ == "__main__":
    if sys.version_info < (3, 6):
        print("错误: 需要 Python 3.6 或更高版本 (支持 hashlib.shake_*)")
    else:
        run_tests()