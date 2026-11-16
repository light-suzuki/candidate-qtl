"""Prepare or complete an immutable A/B review session; no model execution is fabricated."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

from candidate_qtl.claim_review_session import complete_session, finish_wet_review, handoff_sessions, prepare_session


def _json_object(path: str | None) -> dict | None:
    if not path:
        return None
    source = Path(path).resolve()
    raw = source.read_bytes()
    value = json.loads(raw.decode("utf-8-sig"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return {
        **value,
        "artifact_path": str(source),
        "artifact_sha256": hashlib.sha256(raw).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare")
    prepare.add_argument("--route", choices=["a", "b"], required=True)
    prepare.add_argument("--review-policy", choices=["content_v1", "initial_bibliography_v1"], default="content_v1")
    for name in ("hypotheses", "sources", "proposal", "records", "out", "proposer-id", "proposer-model"):
        prepare.add_argument("--" + name, required=True)
    for name in ("bundle", "source-verification", "case-id", "focal-species"):
        prepare.add_argument("--" + name)
    complete = commands.add_parser("complete")
    for name in ("session", "request-sha256"):
        complete.add_argument("--" + name, required=True)
    for name in ("response", "reviewer-id", "reviewer-model"):
        complete.add_argument("--" + name, default="")
    handoff = commands.add_parser("handoff")
    for name in ("a-session", "b-session", "genes", "proteins", "atlas", "out"):
        handoff.add_argument("--" + name, required=True)
    handoff.add_argument("--reference-function")
    handoff.add_argument("--family-index", help="Reference-genome family index TSV for Route A; overrides registry/env/cache.")
    handoff.add_argument("--resources", help="Optional resources.yaml used to resolve an existing family index.")
    handoff.add_argument("--route-c-members")
    handoff.add_argument("--route-c-manifest")
    handoff.add_argument("--identity-assessment", help="Hash-bound post-freeze identity review for derived member Tier evidence; support stays unchanged.")
    handoff.add_argument("--route-b-execution-manifest")
    handoff.add_argument("--route-c-status")
    handoff.add_argument("--allow-unrecorded-route-b", action="store_true")
    handoff.add_argument("--allow-provisional-wet", action="store_true")
    handoff.add_argument("--max-adaptive-rounds", type=int, default=1)
    final = commands.add_parser("finish-wet-review")
    for name in ("handoff", "completion-sha256", "out", "proposals"):
        final.add_argument("--" + name, required=True)
    for name in ("reviewer-id", "reviewer-model"):
        final.add_argument("--" + name, default="")
    final.add_argument("--verification")
    final.add_argument("--independent-review")
    final.add_argument("--executor-id", default="")
    final.add_argument("--executor-model", default="")
    final.add_argument("--allow-provisional-wet", action="store_true")
    final.add_argument("--provisional-screen-ledger")
    final.add_argument("--provisional-screen-ledger-sha256")
    args = parser.parse_args()
    if args.command == "prepare":
        inputs = {k: getattr(args, k) for k in ("hypotheses", "sources", "proposal", "records")}
        for name in ("bundle", "source_verification"):
            if getattr(args, name):
                inputs[name] = getattr(args, name)
        result = prepare_session(args.out, route=args.route, inputs=inputs,
                                 proposer_id=args.proposer_id, proposer_model=args.proposer_model,
                                 case_id=args.case_id or "", focal_species=args.focal_species or "", review_policy=args.review_policy)
    elif args.command == "complete":
        result = complete_session(args.session, args.request_sha256, response_arg=args.response,
                                  reviewer_id=args.reviewer_id, reviewer_model=args.reviewer_model)
    elif args.command == "handoff":
        route_b_execution = _json_object(args.route_b_execution_manifest)
        route_c_execution = _json_object(args.route_c_status)
        if route_c_execution is None and args.route_c_members:
            member_path = Path(args.route_c_members).resolve()
            member_raw = member_path.read_bytes()
            member_count = sum(1 for _ in csv.DictReader(io.StringIO(member_raw.decode("utf-8-sig")), delimiter="\t"))
            route_c_execution = {
                "state": "completed_empty" if member_count == 0 else "completed_members",
                "artifact_path": str(member_path),
                "artifact_sha256": hashlib.sha256(member_raw).hexdigest(),
            }
        result = handoff_sessions(args.a_session, args.b_session, args.out, genes_arg=args.genes,
                                 proteins_arg=args.proteins, atlas_arg=args.atlas, reference_function_arg=args.reference_function,
                                 family_index_arg=args.family_index, resources_arg=args.resources,
                                 route_c_members_arg=args.route_c_members, route_c_manifest_arg=args.route_c_manifest,
                                 identity_assessment_arg=args.identity_assessment,
                                 route_b_execution=route_b_execution, route_c_execution=route_c_execution,
                                 allow_unrecorded_route_b=args.allow_unrecorded_route_b,
                                 allow_provisional_wet=args.allow_provisional_wet,
                                 max_adaptive_rounds=args.max_adaptive_rounds)
    else:
        result = finish_wet_review(args.handoff, args.completion_sha256, args.out, proposals_arg=args.proposals,
                                  reviewer_id=args.reviewer_id, reviewer_model=args.reviewer_model, verification_arg=args.verification,
                                  independent_review_arg=args.independent_review,
                                  executor_id=args.executor_id, executor_model=args.executor_model,
                                  allow_provisional_wet=args.allow_provisional_wet,
                                  provisional_screen_ledger_arg=args.provisional_screen_ledger,
                                  provisional_screen_ledger_sha256=args.provisional_screen_ledger_sha256)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
