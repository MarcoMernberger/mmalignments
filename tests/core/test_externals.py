import subprocess
from mmalignments.models.externals import External  # type: ignore
from mmalignments.models.aligners.bwamem2 import BWAMem2  # type: ignore


def test_build_cmd_and_formatting():
    e = External(
        name="toy",
        primary_binary="toybin",
        parameters={"-a": 1, "-b": True, "--list": ["x", "y"]},
    )

    cmd = e.build_cmd(["positional"], **{"-a": 2})
    # deterministic ordering: --list, -a, -b  (sorted by key strings)
    assert cmd[0] == "toybin"
    assert "-a" in cmd
    assert "2" in cmd
    assert "-b" in cmd
    assert "--list" in cmd


def test_run_with_echo(tmp_path):
    e = External(name="echoer", primary_binary="echo")
    runner = e.run(extra_args=["hello", "world"], capture_output=True)
    # run() now returns a zero-arg callable and executes immediately;
    # the last result is stored on runner.last_result
    assert callable(runner)
    assert hasattr(runner, "last_result")
    cp = runner.last_result
    assert cp is not None
    assert cp.returncode == 0
    assert "hello world" in (cp.stdout or "")

    # calling the returned callable executes again and returns a CompletedProcess
    cp2 = runner()
    assert cp2.returncode == 0


def test_ensure_binary_and_get_version():
    e = External(name="truebin", primary_binary="true")
    assert e.ensure_binary() is True
    # 'true' has no --version; fallback should be returned
    assert e.get_version(fallback="unknown") in ("unknown",)


def test_missing_binary_get_version():
    e = External(name="nope", primary_binary="some_nonexistent_binary_12345")
    assert e.ensure_binary() is False
    assert e.get_version(fallback="fb") == "fb"


def test_bwamem2_align_accepts_post_list_and_runs_in_order(monkeypatch):

    calls = []

    def a():
        calls.append("a")

    def b():
        calls.append("b")

    captured = {}

    def fake_run(*args, **kwargs):
        # align() should wrap the list into a single post-callable
        post_wrapper = kwargs.get("post")
        assert callable(post_wrapper)
        # invoke the wrapper with a fake CompletedProcess
        post_wrapper(subprocess.CompletedProcess(args=["cmd"], returncode=0))
        captured["ran"] = True

        # return a zero-arg callable to mimic External.run behaviour
        return lambda: subprocess.CompletedProcess(args=["cmd"], returncode=0)

    monkeypatch.setattr(BWAMem2, "run", fake_run)

    aln = BWAMem2()
    runner = aln.align(
        index_prefix="/tmp/index", fastq_r1="r1.fq", output_bam="out.bam", post=[a, b]
    )
    assert captured.get("ran") is True
    assert calls == ["a", "b"]


def test_run_writes_header_and_parameters(tmp_path):
    """Verify the combined logfile starts with the RUN LOG header,
    lists run() parameters one-per-line and contains the COMMAND line.
    """
    e = External(name="logger_test", primary_binary="echo")

    # execute with an explicit cwd so the logfile is written into tmp_path
    runner = e.run(
        extra_arguments=["hello", "world"],
        cwd=tmp_path,
        capture_output=True,
        check=True,
        timeout=5,
    )

    # find combined logfile written to cwd
    logs = list(tmp_path.glob("logger_test_*.log"))
    assert logs, "combined logfile not found in cwd"
    content = logs[0].read_text(encoding="utf-8")

    # header and parameter assertions
    assert "RUN LOG" in content.splitlines()[1]
    assert "extra_arguments" in content
    assert "cwd" in content
    assert "capture_output: True" in content
    assert "timeout: 5" in content

    # command logged before subprocess
    assert "COMMAND: echo hello world" in content
