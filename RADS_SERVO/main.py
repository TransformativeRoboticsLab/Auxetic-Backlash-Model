import argparse
import time
from dataclasses import dataclass
from typing import Iterable, List, Optional

try:
    import serial
    from serial.tools import list_ports
except ImportError as exc:
    raise SystemExit(
        "pyserial is required. Install with: pip install pyserial"
    ) from exc


X = 3
Y = 3
EXPECTED_VALUES = 18  # SERVONUM(2) * X(3) * Y(3)


@dataclass
class MatrixPosition:
    linear: List[List[int]]
    expanding: List[List[int]]
    name: str = "unnamed"


def validate_3x3(matrix: List[List[int]], label: str) -> None:
    if len(matrix) != X or any(len(row) != Y for row in matrix):
        raise ValueError(f"{label} must be a {X}x{Y} matrix")


def flatten_3x3(matrix: List[List[int]]) -> List[int]:
    return [int(value) for row in matrix for value in row]


def encode_packet(position: MatrixPosition) -> str:
    validate_3x3(position.linear, "linear")
    validate_3x3(position.expanding, "expanding")

    values = flatten_3x3(position.linear) + flatten_3x3(position.expanding)
    if len(values) != EXPECTED_VALUES:
        raise ValueError(f"Expected {EXPECTED_VALUES} values, got {len(values)}")

    # Packet format: S,<18 values>,E\n
    tokens = ["S", *[str(v) for v in values], "E"]
    return ",".join(tokens) + "\n"


def read_for_window(ser: serial.Serial, seconds: float) -> List[str]:
    lines: List[str] = []
    end_time = time.time() + max(0.0, seconds)
    while time.time() < end_time:
        if ser.in_waiting > 0:
            raw = ser.readline().decode("utf-8", errors="replace").strip()
            if raw:
                lines.append(raw)
        else:
            time.sleep(0.01)
    return lines


def run_sequence(
    ser: serial.Serial,
    positions: Iterable[MatrixPosition],
    dwell_s: float,
    cycles: Optional[int],
    read_window_s: float,
) -> None:
    position_list = list(positions)
    if not position_list:
        raise ValueError("No positions provided")

    cycle_count = 0
    while cycles is None or cycle_count < cycles:
        cycle_count += 1
        print(f"\n=== Cycle {cycle_count} ===")

        for index, position in enumerate(position_list, start=1):
            packet = encode_packet(position)
            ser.write(packet.encode("utf-8"))
            ser.flush()

            print(f"Sent position {index}/{len(position_list)}: {position.name}")
            replies = read_for_window(ser, read_window_s)
            ack_received = False
            nack_messages: List[str] = []
            for line in replies:
                print(f"  Arduino: {line}")
                normalized = line.strip()
                if normalized == "ACK":
                    ack_received = True
                elif normalized.startswith("NACK"):
                    nack_messages.append(normalized)

            if nack_messages:
                print(f"  Result: FAIL ({'; '.join(nack_messages)})")
            elif ack_received:
                print("  Result: OK (ACK)")
            else:
                print("  Result: UNKNOWN (no ACK/NACK seen)")

            time.sleep(max(0.0, dwell_s))


def template_positions() -> List[MatrixPosition]:
    return [
        MatrixPosition(
            name="all_linear_down",
            linear=[
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
            ],
            expanding=[
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
            ],
        ),
        MatrixPosition(
            name="all_linear_up_10",
            linear=[
                [10, 10, 10],
                [10, 10, 10],
                [10, 10, 10],
            ],
            expanding=[
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 0],
            ],
        ),
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Send 3x3 matrix states over serial to RADS_SERVO.ino"
    )
    parser.add_argument(
        "--port",
        default="COM5",
        help="Serial port (default: COM5)",
    )
    parser.add_argument(
        "--baud", type=int, default=9600, help="Baud rate (default: 9600)"
    )
    parser.add_argument(
        "--dwell", type=float, default=1.0, help="Seconds to wait between positions"
    )
    parser.add_argument(
        "--cycles",
        type=int,
        default=1,
        help="How many times to loop positions; use -1 for infinite",
    )
    parser.add_argument(
        "--read-window",
        type=float,
        default=0.35,
        help="Seconds to read Arduino responses after each send",
    )
    parser.add_argument(
        "--list-ports",
        action="store_true",
        help="List available serial ports and exit",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.list_ports:
        ports = list(list_ports.comports())
        if not ports:
            print("No serial ports detected.")
            return
        print("Available ports:")
        for p in ports:
            print(f"  {p.device} - {p.description}")
        return

    cycles: Optional[int] = None if args.cycles < 0 else args.cycles
    positions = template_positions()

    print(f"Opening {args.port} at {args.baud} baud...")
    with serial.Serial(args.port, args.baud, timeout=0.1) as ser:
        time.sleep(2.0)
        startup = read_for_window(ser, 0.5)
        for line in startup:
            print(f"Arduino startup: {line}")

        run_sequence(
            ser=ser,
            positions=positions,
            dwell_s=args.dwell,
            cycles=cycles,
            read_window_s=args.read_window,
        )


if __name__ == "__main__":
    main()
